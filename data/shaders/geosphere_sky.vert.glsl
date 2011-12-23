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

	vec3 eyepos = vec3(varyingEyepos);
	float skyNear, skyFar;
	sphereEntryExitDist(skyNear, skyFar, geosphereCenter, eyepos, geosphereAtmosTopRad);
	vec3 eyedir = normalize(eyepos);
	vec3 a0 = (skyNear * eyedir - geosphereCenter) / geosphereRadius;
	vec3 b0 = (skyFar * eyedir - geosphereCenter) / geosphereRadius;

	vec4 a = vec4(a0.x,a0.y,a0.z,0.0);
	vec4 b = vec4(b0.x,b0.y,b0.z,0.0);

	const float ADF = 500.0;
	float surfaceDensity = atmosColor.w*geosphereAtmosFogDensity;
	vec4 c = vec4(5.8,13.5,33.1,0)*geosphereScale*geosphereRadius/1000000.0;
	//for (int c=0; c<3; c++)
	//ffac[c] = pow(1.0+float(2-c)/6.0,-4);
	const float silly = 2.0*exp(-0.0631-0.418+8.47-15.4+10.5)/ADF;
	const int SN = 8;
	vec4 d = (b-a)/float(SN);
	float len = length(d);

	gl_TexCoord[2] = vec4(0.0,0.0,0.0,1.0);
	mat4 lightDir;
	for (int i=0; i<NUM_LIGHTS; ++i)
		lightDir[i] = normalize(gl_LightSource[i].position - (geosphereCenter.x,geosphereCenter.y,geosphereCenter.z,0.0));

	vec4 occulted = vec4(equal(occultedLight, ivec4(0,1,2,3)));
	vec4 unocculted = vec4(notEqual(occultedLight, ivec4(0,1,2,3)));
	float trad = srad+lrad;
	float absdiff = abs(srad-lrad);

	// estimate integral of scattering along the eyeline
	mat4 scatterInt = mat4(0.0);
	float  scatAtmosInt = 0.0;
	for (int j=0; j<SN+1; j++) {
		vec4 p = a+float(j)*d;
		float r = length(p);
		if (r < 1.0)
			continue;

		float e = exp(-ADF*(r-1.0));

		float lenInvSq = 1.0/dot(p,p);

		vec4 sp = p*(1.0+1.0/ADF);
		float sr = length(sp);
		float se = exp(-ADF*(sr-1.0));

		float slenInvSq = 1.0/dot(sp,sp);

		vec4 y = p * lightDir;
		vec4 a = vec4(1.0) - abs(y)/vec4(r);
		vec4 lt = vec4(lessThan(y,0.0));
		vec4 lightAtmosInt = silly*lt + (1.0-2.0*lt) *
			exp(-0.0631 + a*(-0.418 + a*(8.47 + a*(-15.4 + a*10.5)))) / ADF;
		lightAtmosInt *= vec4(e);

		if (j>0) scatAtmosInt += e * len / 2.0;
		mat4 primaryScatter;
		primaryScatter[0] = exp(-(lightAtmosInt[0]+scatAtmosInt)*c) * e;
		primaryScatter[1] = exp(-(lightAtmosInt[1]+scatAtmosInt)*c) * e;
		primaryScatter[2] = exp(-(lightAtmosInt[2]+scatAtmosInt)*c) * e;
		primaryScatter[3] = exp(-(lightAtmosInt[3]+scatAtmosInt)*c) * e;
		if (j>0) scatAtmosInt += e * len / 2.0;

		vec4 d = y*lenInvSq + sqrt((1.0-lenInvSq)*(1.0-(y*y*lenInvSq)));
		vec4 lightIntensity = clamp(d / (2.0*lightDiscRadii) + 0.5, 0.0, 1.0);
		vec4 projectedLengths;
		if (occultedLight == 0) {
			float dist = length(p - occultCentre - y[0]*lightDir[0] );
			lightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
		}
		else if (occultedLight == 1) {
			float dist = length(p - occultCentre - y[1]*lightDir[0] );
			lightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
		}
		else if (occultedLight == 2) {
			float dist = length(p - occultCentre - y[2]*lightDir[0] );
			lightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
		}
		else if (occultedLight == 3) {
			float dist = length(p - occultCentre - y[3]*lightDir[0] );
			lightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
		}
		primaryScatter[0] *= lightIntensity[0];
		primaryScatter[1] *= lightIntensity[1];
		primaryScatter[2] *= lightIntensity[2];
		primaryScatter[3] *= lightIntensity[3];

		vec4 sy = sp * lightDir;
		vec4 sa = 1.0 - abs(sy)/sr;
		ivec4 slt = lessThan(sy,0.0);
		vec4 secondaryScatterInt = silly*slt + (1.0-2.0*slt) *
			exp(-0.0631 + a*(-0.418 + a*(8.47 + a*(-15.4 + a*10.5)))) / ADF;
		secondaryScatterInt *= se;

		mat4 secondaryScatter;
		secondaryScatter[0] = exp(-(secondaryScatterInt[0]+scatAtmosInt)*c) * se;
		secondaryScatter[1] = exp(-(secondaryScatterInt[1]+scatAtmosInt)*c) * se;
		secondaryScatter[2] = exp(-(secondaryScatterInt[2]+scatAtmosInt)*c) * se;
		secondaryScatter[3] = exp(-(secondaryScatterInt[3]+scatAtmosInt)*c) * se;

		vec4 sd = sy*slenInvSq + sqrt((1.0-slenInvSq)*(1.0-(sy*sy*slenInvSq)));
		vec4 slightIntensity = clamp(sd / (2.0*lightDiscRadii) + 0.5, 0.0, 1.0);
		if (occultedLight == 0) {
			vec3 projectedPoint = sp - sy*lightDir[0];
			float dist = length(projectedPoint - occultCentre);
			slightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ), 0.0, 1.0)));
		}
		else if (occultedLight == 1) {
			vec3 projectedPoint = sp - sy*lightDir[1];
			float dist = length(projectedPoint - occultCentre);
			slightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ), 0.0, 1.0)));
		}
		else if (occultedLight == 2) {
			vec3 projectedPoint = sp - sy*lightDir[2];
			float dist = length(projectedPoint - occultCentre);
			slightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ), 0.0, 1.0)));
		}
		else if (occultedLight == 3) {
			vec3 projectedPoint = sp - sy*lightDir[3];
			float dist = length(projectedPoint - occultCentre);
			slightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ), 0.0, 1.0)));
		}
		secondaryScatter[0] *= slightIntensity[0];
		secondaryScatter[1] *= slightIntensity[1];
		secondaryScatter[2] *= slightIntensity[2];
		secondaryScatter[3] *= slightIntensity[3];

		/*
		   float sy = dot(sp,lightDir[i]);
		   float sa = 1.0 - abs(sy)/sr;
		   float secondaryScatterInt = exp(-0.00287 + sa*(0.459 + sa*(3.83 + sa*(-6.80 + sa*5.25)))) / ADF;
		   if (sy < 0.0) {
		   secondaryScatterInt *= -1;
		   secondaryScatterInt += silly;
		   }
		   secondaryScatterInt *= se;
		   vec3 secondaryScatter = exp(-(secondaryScatterInt+scatAtmosInt[i])*c) * se;

		   float sd = sy*slenInvSq + sqrt((1-slenInvSq)*(1-(sy*sy*slenInvSq)));
		   float slightIntensity = clamp(sd / (2.0*lightDiscRadii[i]) + 0.5, 0.0, 1.0);
		   if (occultedLight == i) {
		   vec3 projectedPoint = sp - sy*lightDir[i];
		   float dist = length(projectedPoint - occultCentre);
		   slightIntensity *= (1.0 - mix(0.0, maxOcclusion,
		   clamp( ( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ), 0.0, 1.0)));
		   }
		   secondaryScatter *= slightIntensity;
		   */

		scatterInt += ( primaryScatter + secondaryScatter ) *
			( (j==0||j==SN) ? 1.0 :
			  2.0 );
	}
	for (int i=0; i<NUM_LIGHTS; ++i) {
		scatterInt[i] *= c*gl_LightSource[i].diffuse*len/2.0;
		for (int c=0; c<3; c++)
			gl_TexCoord[2][c] += scatterInt[i][c];
	}
}
