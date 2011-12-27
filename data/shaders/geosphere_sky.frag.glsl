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

uniform int useSecondary;

varying vec4 varyingEyepos;

#define PI 3.1415926535897931

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

void main(void) {
	vec4 atmosDiffuse = vec4(0.0,0.0,0.0,1.0);
	vec3 eyedir = normalize(vec3(varyingEyepos));
	vec4 mix = 0.0;

	vec4 c = vec4(5.8,13.5,33.1,0)*geosphereScale*geosphereRadius/1000000.0;
	float mc = 2.0*geosphereScale*geosphereRadius/100000.0;
	float mce = mc/0.9;
	float g = 0.76;

	for (int i=0; i<NUM_LIGHTS; ++i) {
		vec3 lightDir = normalize(vec3(gl_LightSource[i].position) - geosphereCenter);
		float mu = dot(eyedir,lightDir);
		mix += gl_TexCoord[3][i] * gl_LightSource[i].diffuse *
		    ( c * (1.0+mu*mu)*(3.0/16*PI) +
		      mc * (3.0/(8.0*PI)) * ( (1-g*g)*(1+mu*mu) ) / ( (2+g*g)*pow(1+g*g-2*g*mu, 1.5)) );

	}
	atmosDiffuse = gl_TexCoord[2] * mix;
	if (useSecondary == 1)
		atmosDiffuse += gl_TexCoord[4];
	atmosDiffuse.a = 1.0;
	//float sun = max(0.0, dot(normalize(eyepos),normalize(vec3(gl_LightSource[0].position))));
	gl_FragColor = atmosDiffuse;

#ifdef ZHACK
	SetFragDepth(gl_TexCoord[6].z);
#endif
}
