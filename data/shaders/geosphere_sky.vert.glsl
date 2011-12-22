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
	vec3 a = (skyNear * eyedir - geosphereCenter) / geosphereRadius;
	vec3 b = (skyFar * eyedir - geosphereCenter) / geosphereRadius;
	vec3 ffac;
	for (int c=0; c<3; c++)
		ffac[c] = pow(1.0+float(2-c)/6.0,-4);
	const float ADF = 500.0;
	float surfaceDensity = atmosColor.w*geosphereAtmosFogDensity;
	gl_TexCoord[2] = vec4(0.0,0.0,0.0,1.0);
	for (int i=0; i<NUM_LIGHTS; ++i) {
		vec3 lightDir = normalize(vec3(gl_LightSource[i].position) - geosphereCenter);

		// estimate integral of scattering along the eyeline
		const int SN = 8;
		vec3 d = (b-a)/float(SN);
		vec3 scatterInt = vec3(0.0,0.0,0.0);
		float atmosDensityTimesRad = 0.0;
		float lastAtmosDensityTimesRad, lightAtmosInt;
		float secondaryScatterInt, secondaryDensity, secondaryScatter;
		float scatAtmosInt = 0.0;
		for (int j=0; j<SN+1; j++) {
			vec3 p = a+j*d;
			float r = length(p);
			if (r < 1.0)
				continue;
			float e = exp(-ADF*(r-1.0));
			lastAtmosDensityTimesRad = atmosDensityTimesRad;
			atmosDensityTimesRad = geosphereScale*geosphereRadius*surfaceDensity*e;
			float y = dot(p,lightDir);
			float a = 1.0 - abs(y)/r;
			lightAtmosInt = e * exp(-0.00287 + a*(0.459 + a*(3.83 + a*(-6.80 + a*5.25)))) / ADF;
			if (y < 0.0)
				lightAtmosInt = 2.0 * exp(-ADF * (sqrt(r*r-y*y)-1.0)) * 0.030854 - lightAtmosInt;
			lightAtmosInt *= geosphereScale*geosphereRadius*surfaceDensity;

			if (j>0)
				scatAtmosInt += (lastAtmosDensityTimesRad + atmosDensityTimesRad) * r/2.0;

			vec3 sp = p*(1.0+2.0/ADF);
			float sr = length(sp);
			float se = exp(-ADF*(sr-1.0));
			float satmosDensityTimesRad = geosphereScale*geosphereRadius*surfaceDensity*se;
			float sy = dot(sp,lightDir);
			float sa = 1.0 - abs(sy)/sr;
			secondaryScatterInt = se * exp(-0.00287 + sa*(0.459 + sa*(3.83 + sa*(-6.80 + sa*5.25)))) / ADF;
			// XXX: commented out just to cut down on instructions!
			//if (sy < 0.0)
				//secondaryScatterInt = 2.0 * exp(-ADF * (sqrt(sr*sr-sy*sy)-1.0)) * 0.030854 - secondaryScatterInt;
			secondaryScatterInt *= geosphereRadius*surfaceDensity;
			secondaryScatter = exp(-(secondaryScatterInt+scatAtmosInt)*ffac) * satmosDensityTimesRad;

			// XXX: this too! I hate shaders!
			//float lightIntensity = clamp(y / (2.0*lightDiscRadii[0]) + 0.5, 0.0, 1.0);

			scatterInt += ( exp(-(lightAtmosInt+scatAtmosInt)*ffac) * atmosDensityTimesRad +
					secondaryScatter) *
				( (j==0||j==SN) ? 1.0 :
				  2.0 );
		}
		scatterInt *= ffac*gl_LightSource[i].diffuse*length(d)/2.0;
		for (int c=0; c<3; c++)
			gl_TexCoord[2][c] += scatterInt[c];
	}
}
