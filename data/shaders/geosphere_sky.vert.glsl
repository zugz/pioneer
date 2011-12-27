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
uniform vec4 occultCentre;
uniform float srad;
uniform float lrad;
uniform float maxOcclusion;

uniform int useSecondary;

varying vec4 varyingEyepos;

varying vec4 rCol;
varying mat4 mCol;
varying vec4 secondaryCol;

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

	// float surfaceDensity = atmosColor.w*geosphereAtmosFogDensity;

	// Numbers for Earth taken from Precomputed Atmospheric Scattering,
	// Bruneton and Neyret 2008
	const float rADF = 795.0;
	const float mADF = 5295.0;

	float scale = geosphereScale*geosphereRadius;

	// Rayleigh scattering for gases, Mie for aerosols.
	// Note: surface density is part of rc and mc; we probably want to
	// separate them out to deal with other planets.
	vec4 rc = scale*vec4(5.8,13.5,33.1,0)/1000000.0;
	vec4 rAbsorption = vec4(0.0);
	vec4 rextinction = rc + rAbsorption;

	float mc = scale*2.0/100000.0;
	float mAbsorption = 0.11*mc;
	float mextinction = mc + mAbsorption;

	// polynomial coefficients calculated by polynomial fitting from ADF
	// values:
	const vec4 rcoeffs = vec4(-0.81786, 10.53578, -19.373, 12.943);
	const float rconstcoeff = 0.10220;
	const vec4 mcoeffs = vec4( -3.808514, 25.688546, -48.801948, 29.360771 );
	const float mconstcoeff = 1.464499;

	// estimated atmosphere density integrals for whole line tangent to
	// sphere:
	const float rTangentInt = 2.0*exp(rconstcoeff+rcoeffs[0]+rcoeffs[1]+rcoeffs[2]+rcoeffs[3])/rADF;
	const float mTangentInt = 2.0*exp(mconstcoeff+mcoeffs[0]+mcoeffs[1]+mcoeffs[2]+mcoeffs[3])/mADF;

	mat4 lightDir = mat4(0.0);
	for (int i=0; i<NUM_LIGHTS; ++i)
		lightDir[i] = normalize(gl_LightSource[i].position - (geosphereCenter.x,geosphereCenter.y,geosphereCenter.z,0.0));

	float trad = srad+lrad;
	float absdiff = abs(srad-lrad);

	// estimate integral of scattering along the eyeline
	mat4 rScatterInt = mat4(0.0);
	mat4 mScatterInt = mat4(0.0);
	vec4 secondaryScatterInt = vec4(0.0);
	float rScatAtmosInt = 0.0;
	float mScatAtmosInt = 0.0;

	const int N = 10;
	vec4 d = (b-a)/float(N);
	float len = length(d);

	for (int j=0; j<N+1; j++) {
		float simpson = (j==0||j==N) ? 1.0 : (j == j/2*2) ? 2.0 : 4.0;
		vec4 p = a+float(j)*d;
		float r = length(p);
		if (r < 1.0)
			continue;

		float re = exp(-rADF*(r-1.0));
		float me = exp(-mADF*(r-1.0));

		float lenInvSq = 1.0/dot(p,p);

		vec4 y = p * lightDir;
		vec4 a = vec4(1.0) - abs(y)/vec4(r);
		vec4 lt = vec4(lessThan(y,vec4(0.0)));
		vec4 rLightAtmosInt =
			re * ( rTangentInt*lt + (1.0-2.0*lt) *
					exp(rconstcoeff + a*(rcoeffs[0] + a*(rcoeffs[1] + a*(rcoeffs[2] + a*rcoeffs[3])))) / rADF);
		vec4 mLightAtmosInt =
			me * ( mTangentInt*lt + (1.0-2.0*lt) *
					exp(mconstcoeff + a*(mcoeffs[0] + a*(mcoeffs[1] + a*(mcoeffs[2] + a*mcoeffs[3])))) / mADF);

		if (j>0) rScatAtmosInt += re * len / 2.0;
		if (j>0) mScatAtmosInt += me * len / 2.0;
		mat4 attenuation;
		// We hand-unroll these loops over lights, to avoid compiler problems:
		attenuation[0] = exp(-( (rLightAtmosInt[0]+rScatAtmosInt)*rextinction + (mLightAtmosInt[0]+mScatAtmosInt)*mextinction ));
		attenuation[1] = exp(-( (rLightAtmosInt[1]+rScatAtmosInt)*rextinction + (mLightAtmosInt[1]+mScatAtmosInt)*mextinction ));
		attenuation[2] = exp(-( (rLightAtmosInt[2]+rScatAtmosInt)*rextinction + (mLightAtmosInt[2]+mScatAtmosInt)*mextinction ));
		attenuation[3] = exp(-( (rLightAtmosInt[3]+rScatAtmosInt)*rextinction + (mLightAtmosInt[3]+mScatAtmosInt)*mextinction ));
		if (j<N) rScatAtmosInt += re * len / 2.0;
		if (j<N) mScatAtmosInt += me * len / 2.0;

		vec4 d = y*lenInvSq + sqrt((1.0-lenInvSq)*(1.0-(y*y*lenInvSq)));
		vec4 lightIntensity = clamp(d / (2.0*lightDiscRadii) + 0.5, 0.0, 1.0);

		// TODO: estimate secondaryLightIntensity by considering sphere around
		// point?
		vec4 secondaryLightIntensity = clamp(d / (6.0*lightDiscRadii) + 0.5, 0.0, 1.0);

		if (occultedLight == 0) {
			float dist = length(p - occultCentre - y[0]*lightDir[0] );
			lightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
			secondaryLightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( 3.0 * ( trad-absdiff ) ), 0.0, 1.0)));
		}
		else if (occultedLight == 1) {
			float dist = length(p - occultCentre - y[1]*lightDir[1] );
			lightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
			secondaryLightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( 3.0 * ( trad-absdiff ) ), 0.0, 1.0)));
		}
		else if (occultedLight == 2) {
			float dist = length(p - occultCentre - y[2]*lightDir[2] );
			lightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
			secondaryLightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( 3.0 * ( trad-absdiff ) ), 0.0, 1.0)));
		}
		else if (occultedLight == 3) {
			float dist = length(p - occultCentre - y[3]*lightDir[3] );
			lightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
			secondaryLightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( 3.0 * ( trad-absdiff ) ), 0.0, 1.0)));
		}

		if (useSecondary == 1) {
			// Complete hack
			vec4 secondaryScatter = vec4(0.0);
			float rse1 = exp(min(0.0,-(rADF*(r-1.0)-1.0)));
			float mse1 = exp(min(0.0,-(mADF*(r-1.0)-1.0)));
			float rse2 = exp(-(mADF*(r-1.0)+1.0));
			float mse2 = exp(-(mADF*(r-1.0)+1.0));
			secondaryScatter += ( rc * re * rse1*exp(-rse1*(1.0/rADF)*rextinction)*rc +
					mc * me * mse1*exp(-mse1*(1.0/mADF)*mextinction)*mc ) * (r-1.0);
			float outer = 1.0 + 6.0/rADF;
			if (r < outer)
				secondaryScatter += ( rc * re * rse2*exp(-rse2*(1.0/rADF)*rextinction)*rc +
						mc * me * mse2*exp(-mse2*(1.0/mADF)*mextinction)*mc ) * (outer-r);

			secondaryScatterInt += simpson * secondaryScatter * secondaryLightIntensity[0] * attenuation[0] * gl_LightSource[0].diffuse;
			secondaryScatterInt += simpson * secondaryScatter * secondaryLightIntensity[1] * attenuation[1] * gl_LightSource[1].diffuse;
			secondaryScatterInt += simpson * secondaryScatter * secondaryLightIntensity[2] * attenuation[2] * gl_LightSource[2].diffuse;
			secondaryScatterInt += simpson * secondaryScatter * secondaryLightIntensity[3] * attenuation[3] * gl_LightSource[3].diffuse;
		}

		attenuation[0] *= lightIntensity[0];
		attenuation[1] *= lightIntensity[1];
		attenuation[2] *= lightIntensity[2];
		attenuation[3] *= lightIntensity[3];

		mat4 rScatter;
		mat4 mScatter;
		rScatter[0] = attenuation[0] * re;
		rScatter[1] = attenuation[1] * re;
		rScatter[2] = attenuation[2] * re;
		rScatter[3] = attenuation[3] * re;
		mScatter[0] = attenuation[0] * me;
		mScatter[1] = attenuation[1] * me;
		mScatter[2] = attenuation[2] * me;
		mScatter[3] = attenuation[3] * me;
		rScatterInt += rScatter * simpson;
		mScatterInt += mScatter * simpson;
	}

	rCol = 0.0;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		// Actual phase function for Rayleigh scattering is
		// (1 + cos(viewAngle) * (3/(16\pi)), but due to multiple scattering,
		// the dependence on angle is barely present in real skies. As part of
		// our implementation of multiple scattering, therefore, we use
		// the constant phase factor:
		const float rPhase = 1.0/(4.0*PI);
		const float hackFactor = 5.0; // <--- XXX HACK XXX
		rCol += rc * hackFactor * gl_LightSource[i].diffuse * rPhase * rScatterInt[i] * len/3.0;

		// Mie scattering, meanwhile, is highly direction-dependent, so we
		// calculate the phase function in the fragment shader.
		mCol[i] = mc * gl_LightSource[i].diffuse * mScatterInt[i] * len/3.0;
	}
	secondaryCol = secondaryScatterInt * len/3.0;
}
