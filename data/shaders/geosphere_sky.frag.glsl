uniform vec4 atmosColor;
// to keep distances sane we do a nearer, smaller scam. this is how many times
// smaller the geosphere has been made
uniform float geosphereScale;
uniform float geosphereAtmosTopRad;
uniform vec3 geosphereCenter;
uniform float geosphereAtmosFogDensity;

uniform float geosphereRadius;

uniform int useSecondary;

varying vec4 varyingEyepos;

varying vec4 rCol;
varying mat4 mCol;
varying vec4 secondaryCol;

#define PI 3.1415926535897931

void main(void) {
	vec4 atmosDiffuse = vec4(0.0,0.0,0.0,1.0);
	vec3 eyedir = normalize(vec3(varyingEyepos));

	// For the Earth's atmosphere:
	const float g = 0.76;

	atmosDiffuse = rCol;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		vec3 lightDir = normalize(vec3(gl_LightSource[i].position) - geosphereCenter);
		float mu = dot(eyedir,lightDir);
		// taken from Bruneton (loc. cit.):
		float mPhase = (3.0/(8.0*PI)) * ( (1-g*g)*(1+mu*mu) ) / ( (2+g*g)*pow(1+g*g-2*g*mu, 1.5));

		atmosDiffuse += mCol[i] * mPhase;
	}
	if (useSecondary == 1)
		atmosDiffuse += secondaryCol;
	atmosDiffuse.a = 1.0;
	//float sun = max(0.0, dot(normalize(eyepos),normalize(vec3(gl_LightSource[0].position))));
	gl_FragColor = atmosDiffuse;

#ifdef ZHACK
	SetFragDepth(gl_TexCoord[0].z);
#endif
}
