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
	for (int c=0; c<3; c++)
		ffac[c] = pow(1.0+float(2-c)/6.0,-4);
	const float ADF = 500.0;
	float surfaceDensity = atmosColor.w*geosphereAtmosFogDensity;
	gl_TexCoord[2] = vec4(0.0,0.0,0.0,1.0);
	for (int i=0; i<NUM_LIGHTS; ++i) {
		vec3 lightDir = normalize(vec3(gl_LightSource[i].position) - geosphereCenter);

		// estimate integral of scattering along the eyeline
		const int SN = 6;
		vec3 d = (b-a)/float(SN);
		vec3 scatterInt = vec3(0.0,0.0,0.0);
		float atmosDensity, lastAtmosDensity, lightAtmosInt;
		float scatAtmosInt = 0.0;
		for (int j=0; j<SN+1; j++) {
			vec3 p = a+j*d;
			{
				// estimate integral of density
				const int lSN=3;
				float lb = -dot(p,lightDir);
				vec3 ld = ((lb + sqrt( lb*lb - dot(p,p) + 1 ))/float(lSN))*lightDir;
				lightAtmosInt = surfaceDensity*exp(-ADF*(length(p)-0.985));
				for (int k=1; k<lSN; k++) {
					vec3 lp = p + k*ld;
					// 0.985 == 1/ATMOSPHERE_RADIUS
					lightAtmosInt += 2.0*surfaceDensity*exp(-ADF*(length(lp)-0.985));
				}
				lightAtmosInt *= length(ld)*geosphereAtmosTopRad/2.0;
			}

			lastAtmosDensity = atmosDensity;
			atmosDensity = surfaceDensity*exp(-ADF*(length(p)-0.985));

			if (j>0)
				scatAtmosInt += (lastAtmosDensity + atmosDensity) *
					length(d)*geosphereAtmosTopRad/2.0;

			scatterInt += exp(-(lightAtmosInt+scatAtmosInt)*ffac) *
				atmosDensity *
				( (j==0||j==SN) ? 1.0 :
				  (j==1||j==3||j==5) ? 4.0 :
				  2.0 );
		}
		//float lightIntensity = intensityOfLightAtPoint(a+(d*(SN/2))*geosphereRadius/geosphereAtmosTopRad,
				//lightDir, lightDiscRadii[i], occultedLight == i, occultCentre, srad,
				//lrad, maxOcclusion);
		scatterInt *= ffac*gl_LightSource[i].diffuse*length(d)*geosphereAtmosTopRad/3.0;
		for (int c=0; c<3; c++)
			gl_TexCoord[2][c] += scatterInt[c];
	}
}
