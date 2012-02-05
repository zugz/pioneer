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
varying vec4 extraIn;

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


// FIXME: check what happens when geosphereScale > 1

void main(void)
{
#ifdef ZHACK
	gl_Position = logarithmicTransform();
#else
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
#endif
	varyingEyepos = gl_ModelViewMatrix * gl_Vertex;

	vec3 eyepos = vec3(varyingEyepos);

	// float surfaceDensity = atmosColor.w*geosphereAtmosFogDensity;

	float scale = geosphereScale*geosphereRadius;

	// Numbers for Earth taken from Precomputed Atmospheric Scattering,
	// Bruneton and Neyret 2008
	//
	// Rayleigh scattering for gases, Mie for aerosols.
	const vec4 gasScatteringSeaLevel = vec4(5.8e-6,13.5e-6,33.1e-6,0);
	// rADF = planetRadius / gasScaleHeight
	const float rADF = 795.0;
	vec4 rc = scale*gasScatteringSeaLevel;
	vec4 rAbsorption = vec4(0.0);
	vec4 rextinction = rc + rAbsorption;

	const float particleScatteringSeaLevel = 20e-6;
	// rADF = planetRadius / particleScaleHeight
	const float mADF = 5295.0;
	float mc = scale*particleScatteringSeaLevel;
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

	const float constantPhase = 1.0/(4.0*PI);

	mat4 lightDir = mat4(0.0);
	mat4 lightDiffuse = mat4(0.0);
	for (int i=0; i<NUM_LIGHTS; ++i) {
		lightDir[i] = normalize(gl_LightSource[i].position);
		lightDiffuse[i] = gl_LightSource[i].diffuse;
	}

	float trad = srad+lrad;
	float absdiff = abs(srad-lrad);

	// estimate integral of scattering along the line from the eye to the
	// vertex
	mat4 rScatterInt = mat4(0.0);
	mat4 mScatterInt = mat4(0.0);
	extraIn = vec4(0.0);
	float rScatAtmosInt = 0.0;
	float mScatAtmosInt = 0.0;

	vec4 a,b;
	{
		float atmosNear;
		float atmosFar;
		sphereEntryExitDist(atmosNear, atmosFar, geosphereCenter, eyepos, geosphereAtmosTopRad);

		vec3 eyedir = normalize(eyepos);
		a = vec4(atmosNear*eyedir - geosphereCenter, 0.0) / geosphereRadius;
#ifdef GROUND
		b = vec4(eyepos - geosphereCenter, 0.0) / geosphereRadius;
#else
		b = vec4(atmosFar*eyedir - geosphereCenter, 0.0) / geosphereRadius;
#endif
	}

	// N: steps in our numerical integration. Must be even, since we use
	// Simpson's rule. TODO: vary according to distance relative to 1/ADF?
	//
	// More importantly: we can't do (all) ground here, there are just too
	// many damned vertices for something this expensive. I'm currently
	// thinking that we should instead do all ground on the frag shader,
	// though still using vert for sky, with whatever crappy ad hoccities we
	// can find to get aerial perspective working on the cheap while also
	// having approximately correct continuation of atmospheric effects. Then
	// we can use the same tricks for models.
#ifdef GROUND
	const int N = 6;
#else
	const int N = 10;
#endif
	vec4 d = (b-a)/float(N);
	float len = length(d);
	vec4 p = a;
	for (int j=0; j<N+1; j++) {
		float simpson = (j==0||j==N) ? 1.0 : (j == j/2*2) ? 2.0 : 4.0;
		p += d;
		float r = length(p);
#ifndef GROUND
		if (r < 1.0)
			continue;
#endif

		float re = exp(-rADF*(r-1.0));
		float me = exp(-mADF*(r-1.0));

		float lenInvSq = 1.0/(r*r);

		vec4 y = p * lightDir;
		vec4 z = vec4(1.0) - abs(y)/vec4(r);
		vec4 lt = vec4(lessThan(y,vec4(0.0)));
		vec4 rLightAtmosInt =
			re * ( rTangentInt*lt + (1.0-2.0*lt) *
					exp(rconstcoeff + z*(rcoeffs[0] + z*(rcoeffs[1] + z*(rcoeffs[2] + z*rcoeffs[3])))) / rADF);
		vec4 mLightAtmosInt =
			me * ( mTangentInt*lt + (1.0-2.0*lt) *
					exp(mconstcoeff + z*(mcoeffs[0] + z*(mcoeffs[1] + z*(mcoeffs[2] + z*mcoeffs[3])))) / mADF);

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

		vec4 ld = y*lenInvSq + sqrt(max(vec4(0.00001),(1.0-lenInvSq)*(vec4(1.0)-(y*y*lenInvSq))));
		vec4 lightIntensity = clamp(ld / (2.0*lightDiscRadii) + 0.5, vec4(0.0), vec4(1.0));

		vec4 secondaryLightIntensity = vec4(1.0);
		const float secMaxDist = 5.0/rADF;
		//secondaryLightIntensity = clamp(ld / (2.0*(lightDiscRadii+secMaxDist)) + 0.5, 0.0, 1.0);

		/*
		if (occultedLight == 0) {
			float dist = length(p - occultCentre - y[0]*lightDir[0] );
			lightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
				secondaryLightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
							clamp( ( trad-dist + secMaxDist ) / ( trad-absdiff + 2.0*secMaxDist ), 0.0, 1.0)));
		}
		else if (occultedLight == 1) {
			float dist = length(p - occultCentre - y[1]*lightDir[1] );
			lightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
				secondaryLightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
							clamp( ( trad-dist + secMaxDist ) / ( trad-absdiff + 2.0*secMaxDist ), 0.0, 1.0)));
		}
		else if (occultedLight == 2) {
			float dist = length(p - occultCentre - y[2]*lightDir[2] );
			lightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
				secondaryLightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
							clamp( ( trad-dist + secMaxDist ) / ( trad-absdiff + 2.0*secMaxDist ), 0.0, 1.0)));
		}
		else if (occultedLight == 3) {
			float dist = length(p - occultCentre - y[3]*lightDir[3] );
			lightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
				secondaryLightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
							clamp( ( trad-dist + secMaxDist ) / ( trad-absdiff + 2.0*secMaxDist ), 0.0, 1.0)));
		}
		*/

		mat4 rScatter;
		mat4 mScatter;
		rScatter[0] = attenuation[0] * lightIntensity[0] * re;
		rScatter[1] = attenuation[1] * lightIntensity[1] * re;
		rScatter[2] = attenuation[2] * lightIntensity[2] * re;
		rScatter[3] = attenuation[3] * lightIntensity[3] * re;
		mScatter[0] = attenuation[0] * lightIntensity[0] * me;
		mScatter[1] = attenuation[1] * lightIntensity[1] * me;
		mScatter[2] = attenuation[2] * lightIntensity[2] * me;
		mScatter[3] = attenuation[3] * lightIntensity[3] * me;
		rScatterInt += rScatter * simpson;
		mScatterInt += mScatter * simpson;

		vec4 secondaryScatter = vec4(0.0);
		if (useSecondary == 1)
		{
		//	extraIn += simpson * (len/3.0) * re * rc * PI * (1.0/(4.0*PI)) * exp(-(rScatAtmosInt*rextinction + mScatAtmosInt*mextinction)) * (lightDiffuse * secondaryLightIntensity);
			//secondaryScatter = rc * re * (0.1/rADF) * 4*PI * (1.0 + rc * re * (0.1/rADF) * 4*PI * (1.0 + rc * re * (0.1/rADF) * 4*PI));

			// Idea: 1.0/rADF is gaseous optical depth for a radial line from
			// the surface to infinity. Our fairly arbitrary estimate is that
			// through every point passes multiply scattered light
			// corresponding to twice the amount scattered along this optical
			// depth, with scattering path being of optical depth that of
			// direct light to the point.
			secondaryScatter = 4.0*PI*rc * 2.0/rADF * exp(-rextinction*1.0/rADF);

			// Consider the cone formed from line segments starting at p,
			// passing tangent to the planet, and ending at the "edge of the
			// atmosphere" (which we arbitrarily consider to be at radius
			// 1+3/rADF). The proportion of this cone which is unshadowed
			// is a decent estimate for the secondary scatter p is receiving,
			// as a proportion of that it does in full uneclipsed daylight.
			//
			// For shadowing by the planet, we estimate the unshadowed
			// proportion of the cone by the unshadowed proportion of its
			// bounding circle:

#ifdef GROUND
			r *= 1.0+1.0/rADF;
			y *= 1.0+1.0/rADF;
#endif
			const float atmoRad = 1.0 + 3.0/rADF;
			float r0;
			if (r < 1.00001) r0 = 1.0;
			else r0 = r - (r-1.0/r)*(sqrt((atmoRad*atmoRad-1.0)/(r*r-1.0))+1.0);
			float atmoHorizonRad = sqrt(atmoRad*atmoRad - r0*r0);

			vec4 sinphi = y/r;
			vec4 cosphi = sqrt(1.0-sinphi*sinphi);
			vec4 s0 = cosphi*r0;
			vec4 posUnshadowed = vec4(greaterThan(sinphi*r0 - atmoHorizonRad*cosphi, vec4(0.0)));
			vec4 negUnshadowed = vec4(greaterThan(sinphi*r0 + atmoHorizonRad*cosphi, vec4(0.0)));

			// quadratic in cos(\theta) - in shadow cylinder of the planet iff
			// (atmoHorizonRad*sin(\theta))^2 + (s0 + atmoHorizonRad*cos(\theta)*sin(\phi))^2 <= 1
			vec4 A2 = max(0.0001, 2.0*vec4(atmoHorizonRad*atmoHorizonRad)*(vec4(1.0)-sinphi*sinphi));
			vec4 B = -2.0*vec4(atmoHorizonRad)*s0*sinphi;
			vec4 C2 = 2.0*(vec4(1.0) - vec4(atmoHorizonRad*atmoHorizonRad)+s0*s0);

			vec4 rootdet = sqrt(max(vec4(0.0), B*B - A2*C2));

			vec4 theta1 = acos(mix(clamp((-B+rootdet)/A2,-1.0,1.0), 1.0, posUnshadowed));
			vec4 theta2 = acos(mix(clamp((-B-rootdet)/A2,-1.0,1.0), -1.0, negUnshadowed));
			secondaryLightIntensity *= (theta2-theta1)/PI;

			/*
			// XXX: Entirely ad hoc secondary scatter: split into two components, from below and
			// from above, and in each case just integrate along a line, using density 1/ADF away as
			// estimate of average... there's no good theoretical reasoning behind this, but it
			// seems to work!
			vec4 secondaryScatter = vec4(0.0);
#ifndef GROUND
			float samp1 = max(0.0,r-1.0/rADF-1.0);
			float rse1 = exp(-(rADF*samp1));
			float mse1 = exp(-(mADF*samp1));
			secondaryScatter += exp(-(rse1*(r-1)*rextinction + mse1*(r-1)*mextinction)) * (rc * re + mc * me) * (r-1.0);
#endif
			float samp2 = r+1.0/rADF-1.0;
			float rse2 = exp(-(rADF*samp2));
			float mse2 = exp(-(mADF*samp2));
			float outer = 1.0 + 6.0/rADF;
			if (r < outer)
				secondaryScatter += exp(-(rse2*(outer-r)*rextinction + mse2*(r-1)*mextinction)) * (rc * re + mc * me) * (outer-r);
			*/

			extraIn += simpson * (len/3.0) * (re * rc + me * mc) * PI * constantPhase * secondaryScatter * (matrixCompMult(attenuation, lightDiffuse) * secondaryLightIntensity);
			// extraIn += simpson * (len/3.0) * (re * rc + me * mc) * PI * constantPhase * secondaryScatter * secondaryLightIntensity[0];
		}

#ifdef GROUND
		if (j==N) {
			vec4 direct = lightIntensity * max(vec4(0.0), vec4(normalize(gl_NormalMatrix * gl_Normal),0.0) * lightDir);
			//extraIn += gl_Color * (direct[0] + secondaryLightIntensity[0] * PI * (1.0/(2.0*PI)) * secondaryScatter) * lightDiffuse[0];
			extraIn += gl_Color * (direct[0] + secondaryLightIntensity[0] * PI * (1.0/(2.0*PI)) * secondaryScatter) * attenuation[0] * lightDiffuse[0];
			extraIn += gl_Color * (direct[1] + secondaryLightIntensity[1] * PI * (1.0/(2.0*PI)) * secondaryScatter) * attenuation[1] * lightDiffuse[1];
			extraIn += gl_Color * (direct[2] + secondaryLightIntensity[2] * PI * (1.0/(2.0*PI)) * secondaryScatter) * attenuation[2] * lightDiffuse[2];
			extraIn += gl_Color * (direct[3] + secondaryLightIntensity[3] * PI * (1.0/(2.0*PI)) * secondaryScatter) * attenuation[3] * lightDiffuse[3];
		}
#endif
	}

	rCol = 0.0;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		// Actual phase function for Rayleigh scattering is (1 + cos^2(viewAngle)) * (3/(16\pi), but
		// due to multiple scattering, the dependence on angle is barely present in real skies. As
		// part of our implementation of multiple scattering, therefore, we use the constant phase
		// factor:
		// float mu = dot(varyingEyepos, lightDir[i]);
		// float rPhase = (1.0+mu*mu) * (3.0/(16.0*PI));
		const float rPhase = 1.0/(4.0*PI);
		// We also just multiply up the value by an ad hoc factor. We could try to excuse this by
		// saying it accounts for multiple scatter, but I've no reason other than subjective
		// judgement of results to think it physically realistic!
		const float rhackFactor = 1.0; // <--- XXX HACK XXX
		rCol += rc * rhackFactor * PI * gl_LightSource[i].diffuse * rPhase * rScatterInt[i] * len/3.0;

		// Mie scattering, meanwhile, is highly direction-dependent, so we
		// calculate the phase function in the fragment shader.
		const float mhackFactor = 1.0; // <--- XXX HACK XXX
		mCol[i] = mc * mhackFactor * PI * gl_LightSource[i].diffuse * mScatterInt[i] * len/3.0;
	}
}
