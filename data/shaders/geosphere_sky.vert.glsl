uniform vec4 atmosColor;
// to keep distances sane we do a nearer, smaller scam. this is how many times
// smaller the geosphere has been made
uniform float geosphereScale;
uniform float geosphereAtmosTopRad;
uniform vec3 geosphereCenter;
uniform float geosphereAtmosFogDensity;

uniform float geosphereRadius;
uniform vec4 lightDiscRadii;
uniform int occultedLight;
uniform vec3 occultCentre;
uniform float srad;
uniform float lrad;
uniform float maxOcclusion;

varying vec4 varyingEyepos;

void sphereEntryExitDist(out float near, out float far, in vec3 sphereCenter, in vec3 eyeTo, in float radius)
{
	vec3 v = -sphereCenter;
	vec3 dir = normalize(eyeTo);
	float b = -dot(v, dir);
	float det = (b * b) - dot(v, v) + (radius * radius);
	near = 0.0;
	far = 0.0;
	if (det > 0.0) {
		det = sqrt(det);
		float i1 = b - det;
		float i2 = b + det;
		if (i2 > 0.0) {
			near = max(i1, 0.0);
			far = i2;
		}
	}
}

void main(void)
{
#ifdef ZHACK
	gl_Position = logarithmicTransform();
#else
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
#endif
	varyingEyepos = gl_ModelViewMatrix * gl_Vertex;

	float skyNear, skyFar;
	sphereEntryExitDist(skyNear, skyFar, geosphereCenter, varyingEyepos, geosphereAtmosTopRad);
	vec3 eyedir = normalize(varyingEyepos);
	vec3 a = (skyNear * eyedir - geosphereCenter) / geosphereAtmosTopRad;
	vec3 b = (skyFar * eyedir - geosphereCenter) / geosphereAtmosTopRad;
	vec3 ffac;
	for (int c=0; c<3; c++) {
		float freq = 1.0/(1.0+float(2-c)/6.0);
		ffac[c] = freq*freq*freq*freq;
	}
	const float ADF = 500.0;
	float surfaceDensity = atmosColor.w*geosphereAtmosFogDensity;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		vec3 lightDir = normalize(vec3(gl_LightSource[i].position) - geosphereCenter);

		// estimate integral of scattering along the eyeline
		const int SN = 2;
		vec3 d = (b-a)/float(SN);
		float atmosSamples[SN+1];
		float lightIntensity[SN+1];
		for (int j=0; j<SN+1; j++) {
			vec3 p = a+j*d;
			atmosSamples[j] = surfaceDensity*exp(-ADF*(length(p)-0.985));
			lightIntensity[j] = intensityOfLightAtPoint(p*geosphereRadius/geosphereAtmosTopRad,
					lightDir, lightDiscRadii[i], occultedLight == i, occultCentre, srad,
					lrad, maxOcclusion);
		}

		vec3 scatterInt = vec3(0.0,0.0,0.0);
		for (int j=0; j<SN+1; j++) {
			vec3 p = a+j*d;
			float lightAtmosInt = 0.0;
			{
				// estimate integral of density
				const int lSN=2;
				float lb = -dot(p,lightDir);
				vec3 ld = ((lb + sqrt( lb*lb - dot(p,p) + 1 ))/float(lSN))*lightDir;
				for (int k=0; k<lSN+1; k++) {
					vec3 lp = p + k*ld;
					// 0.985 == 1/ATMOSPHERE_RADIUS
					lightAtmosInt += surfaceDensity*exp(-ADF*(length(lp)-0.985)) *
						( (k==0||k==lSN) ? 1.0 :
						  (k==1||k==3||k==5) ? 4.0 :
						  2.0 );
				}
				lightAtmosInt *= length(ld)*geosphereAtmosTopRad/3.0;
			}

			float scatAtmosInt = 0.0;
			if (j>0) {
				scatAtmosInt += atmosSamples[0];
				for (int k=1; k<j; k++)
					scatAtmosInt += atmosSamples[k] * 2.0;
				scatAtmosInt += atmosSamples[j];
				scatAtmosInt *= length(d)*geosphereAtmosTopRad/2.0;
			}

			for (int c=0; c<3; c++) {
				float lightFogFactor = exp(-lightAtmosInt*ffac[c]);
				float scatFogFactor = exp(-scatAtmosInt*ffac[c]);
				scatterInt[c] += lightFogFactor * scatFogFactor * ffac[c] * gl_LightSource[i].diffuse[c] *
					atmosSamples[j] * lightIntensity[j] *
					( (j==0||j==SN) ? 1.0 :
					  (j==1||j==3||j==5) ? 4.0 :
					  2.0 );
			}
		}
		for (int c=0; c<3; c++) {
			scatterInt[c] *= length(d)*geosphereAtmosTopRad/3.0;
			gl_TexCoord[2][c] = scatterInt[c];
		}
	}
}
