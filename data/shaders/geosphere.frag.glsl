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
	vec3 eyepos = vec3(gl_TexCoord[1]);
	vec3 tnorm = normalize(vec3(gl_TexCoord[2]));
	vec4 total = vec4(0.0);

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

	const float g = 0.76;

	const float constantPhase = 1.0/(4.0*PI);

	// estimated atmosphere density integrals for whole line tangent to
	// sphere:
	const float rTangentInt = 2.0*exp(rconstcoeff+rcoeffs[0]+rcoeffs[1]+rcoeffs[2]+rcoeffs[3])/rADF;
	const float mTangentInt = 2.0*exp(mconstcoeff+mcoeffs[0]+mcoeffs[1]+mcoeffs[2]+mcoeffs[3])/mADF;

	mat4 lightDir = mat4(0.0);
	mat4 lightDiffuse = mat4(0.0);
	for (int i=0; i<NUM_LIGHTS; ++i) {
		lightDir[i] = normalize(gl_LightSource[i].position);
		lightDiffuse[i] = gl_LightSource[i].diffuse;
	}

	float trad = srad+lrad;
	float absdiff = abs(srad-lrad);

	// when does the eye ray intersect atmosphere
	float atmosDist = findSphereEyeRayEntryDistance(geosphereCenter, eyepos, geosphereAtmosTopRad);

	vec3 eyedir = normalize(eyepos);
	vec4 a = vec4(atmosDist*eyedir - geosphereCenter,0.0) / geosphereRadius;
	vec4 b = vec4(eyepos - geosphereCenter,0.0) / geosphereRadius;

	float r1 = length(a);
	float r2 = length(b);
	float re1 = exp(-rADF*(r1-1.0));
	float me1 = exp(-mADF*(r1-1.0));
	float re2 = exp(-rADF*(r2-1.0));
	float me2 = exp(-mADF*(r2-1.0));

	vec4 ly = b * lightDir;
	vec4 lz = vec4(1.0) - abs(ly)/vec4(r2);
	vec4 llt = vec4(lessThan(ly,vec4(0.0)));
	vec4 rs = ( rTangentInt*llt + (1.0-2.0*llt) *
				exp(rconstcoeff + lz*(rcoeffs[0] + lz*(rcoeffs[1] + lz*(rcoeffs[2] + lz*rcoeffs[3])))) / rADF);
	vec4 ms = ( mTangentInt*llt + (1.0-2.0*llt) *
				exp(mconstcoeff + lz*(mcoeffs[0] + lz*(mcoeffs[1] + lz*(mcoeffs[2] + lz*mcoeffs[3])))) / mADF);
	vec4 rLightAtmosInt = re2 * rs;
	vec4 mLightAtmosInt = me2 * ms;

	float rViewAtmosInt, mViewAtmosInt;
	{
		float vy1 = dot(vec3(a),eyedir);
		float vz1 = 1.0 - abs(vy1)/r1;
		float vy2 = dot(vec3(b),eyedir);
		float vz2 = 1.0 - abs(vy2)/r2;
		// TODO: could probably optimise this when the view point is low
		// enough by using flat earth hypothesis, and corresponding analytic
		// solution to the integral, and when we're in space by noting that
		// one term is 0.
		float vrs1 = exp(rconstcoeff + vz1*(rcoeffs[0] + vz1*(rcoeffs[1] + vz1*(rcoeffs[2] + vz1*rcoeffs[3])))) / rADF;
		float vrs2 = exp(rconstcoeff + vz2*(rcoeffs[0] + vz2*(rcoeffs[1] + vz2*(rcoeffs[2] + vz2*rcoeffs[3])))) / rADF;
		float vms1 = exp(mconstcoeff + vz1*(mcoeffs[0] + vz1*(mcoeffs[1] + vz1*(mcoeffs[2] + vz1*mcoeffs[3])))) / mADF;
		float vms2 = exp(mconstcoeff + vz2*(mcoeffs[0] + vz2*(mcoeffs[1] + vz2*(mcoeffs[2] + vz2*mcoeffs[3])))) / mADF;
		rViewAtmosInt = ((vy1>0.0) == (vy2>0.0)) ? abs(re1*vrs1-re2*vrs2) : (re1+re2)/2.0 * (rTangentInt - (vrs1+vrs2));
		mViewAtmosInt = ((vy1>0.0) == (vy2>0.0)) ? abs(me1*vms1-me2*vms2) : (me1+me2)/2.0 * (mTangentInt - (vms1+vms2));
	}

	// d[i] = dot(lightDir[i],t[i]) where t[i] is the unique point on the unit
	// sphere whose tangent plane contains b, is in the plane of lightDir[i] and
	// b, and is towards the light from b.
	float lenInvSq = 1.0/dot(b,b);
	vec4 d = ly*lenInvSq + sqrt((1.0-lenInvSq)*(1.0-(ly*ly*lenInvSq)));
	vec4 lightIntensity = clamp(d / (2.0*lightDiscRadii) + 0.5, 0.0, 1.0);

	vec4 secondaryScatter = vec4(0.0);
	if (useSecondary == 1)
	{
		secondaryScatter = 4.0*PI * rc * 2.0/rADF * exp(-rextinction*1.0/rADF);
		/*
		float samp1 = max(0.0,r2-1.0/rADF-1.0);
		float samp2 = r2+1.0/rADF-1.0;
		float rse1 = exp(-(rADF*samp1));
		float rse2 = exp(-(rADF*samp2));
		float mse1 = exp(-(mADF*samp1));
		float mse2 = exp(-(mADF*samp2));
		secondaryScatter += exp(-(rse1*(r2-1)*rextinction + mse1*(r2-1)*mextinction)) * (rc * re2 + mc * me2) * (r2-1.0);
		float outer = 1.0 + 6.0/rADF;
		secondaryScatter += exp(-(rse1*(outer-r2)*rextinction + mse2*(r2-1.0)*mextinction)) * (rc * re2 + mc * me2) * (outer-r2);
		*/

		//secondaryScatter = rc * re2 * (0.1/rADF) * 4*PI * (1.0 + rc * re2 * (0.1/rADF) * 4*PI * (1.0 + rc * re2 * (0.1/rADF) * 4*PI));
	}

	//for (int i=0; i<NUM_LIGHTS; ++i) {
	int i=0;
		if (occultedLight == i) {
			float dist = length(b - occultCentre - ly[i]*lightDir[i] );
			lightIntensity[i] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
		}
		vec4 inAttenuation = exp(-(rLightAtmosInt[i]*rextinction + mLightAtmosInt[i]*mextinction));
		vec4 outAttenuation = exp(-(rViewAtmosInt*rextinction + mViewAtmosInt*mextinction));

		float nDotVP = max(0.0, dot(vec4(tnorm,0.0), lightDir[i]));
		total += lightIntensity[i] * inAttenuation * gl_Color * lightDiffuse[i] * (nDotVP + (1.0/(2*PI)) * secondaryScatter) * outAttenuation;

			// aerial perspective
			// We use flat earth hypothesis to make the integral simple.
			// Problem: the numerical integration is *still* too expensive.
			// We could separate the Mie from the Rayleigh, and separately do
			// each integral analytically. This would be quite wrong, but
			// would at least be fairly cheap. Could we can hack in some
			// factor to one to account for the other?
			// Problem even with that: we need secondary scatter to be
			// (affine) linear in exp(-rADF*h) to be able to do the integral.
			// Really need to find a decent such hack for skylight!
			//
			// Alternative to separating out the two particles: see
			// Preetham-Shirley-Smits, where they use a cubic approximation to
			// one of the exponentials hence allowing a by-parts analytic
			// integral. They leave as a factor to be (pre)computed, however,
			// the light scattered at a point at sea level in a particular
			// direction - making the assumption that this scales with density
			// as we go up. That's a reasonable assumption when we're viewing
			// from the ground, less so when we're viewing from space, but
			// perhaps we can get away with it.
			//
			// So what we need, for this and for the numerical integration in
			// the sky shader, is a decent estimate for the secondary scatter
			// at a given point, depending on height and sun-angle. If we
			// follow PSS and use it just as a multiplicative factor here, we
			// don't need to worry about it being a nice integrable function.
			//
			// Idea for a hacky solution to that: calculate proportion of
			// sphere of some arbitrary radius like 1+rADF/5 which is visible
			// from this position and lit, and pretend that the entire radial
			// line passing through each such point is visible and lit;
			// furthermore, pretend for attenuation-calculation purposes that
			// the scattered light all comes to us from the point on this
			// sphere radially above/below us (a systematic underestimate, but
			// probably not significant). Hmm, no, this seems not to lead to
			// pleasant calculations.
			//
			// Am I overcomplicating this? The intensity of secondary-scatter
			// seems like it should be more-or-less constant with height, so
			// perhaps rhackFactor is actually the right solution after all?
			// Just need to have it and the phase function vary reasonably
			// with wavelength, and multiply by the proportion of visible sky
			// which is directly lit. Of course, rhackFactor isn't literally
			// right - we do want the ground to get some secondary lighting
			// too. So how about: main scatter uses rayleigh phase function,
			// then we have a *constant* secondary lighting intensity with
			// constant phase.
			//
			//
			// Another issue: the flat-earth hypothesis does lead to
			// noticeable under-attenuation when we view from space. One hack
			// to try: scale by the ratio between the actual optical depth and
			// the flat-earth optical depth. Alternative if we can get the
			// branching working right: use the sky shader for the horizon
			// from space.
			float D = distance(a,b);
			vec4 rInScatter = vec4(0.0);
			vec4 mInScatter = vec4(0.0);
			if (abs(r2-r1) < 0.05/rADF) {
				vec4 sumext = rextinction*re1 + mextinction*me1;
				//vec4 I = exp(-(rext*rs[i] + mext*ms[i])) * (1.0 - exp(-(rext+mext)*D)) / (rext+mext);
				vec4 I = inAttenuation * (1.0 - exp(-sumext*D)) / sumext;
				rInScatter = rc * re1 * I;
				mInScatter = mc * me1 * I;
			} else {
				const int N=6;
				float dsbydh = D/(r2-r1);
				float Dh = (r2-r1)/float(N);
				float rstep = exp(-rADF*Dh);
				float mstep = exp(-mADF*Dh);
				float reh = re1;
				float meh = me1;
				vec4 rfac = rextinction * (rs[i] - dsbydh/rADF);
				float mfac = mextinction * (ms[i] - dsbydh/mADF);
				for (int j=0; j<N+1; j++) {
					float simpson = (j==0||j==N) ? 1.0 : (j == j/2*2) ? 2.0 : 4.0;
					vec4 atten = exp(-(reh * rfac + meh * mfac));
					rInScatter += simpson * atten * rc*reh;
					mInScatter += simpson * atten * mc*meh;
					reh *= rstep;
					meh *= mstep;
				}
				vec4 atten0 = exp(-dsbydh * (re1 * rextinction / rADF + me1 * mextinction / mADF));
				rInScatter *= dsbydh * atten0 * Dh/3.0;
				mInScatter *= dsbydh * atten0 * Dh/3.0;
			}

			float mu = dot(eyedir, vec3(lightDir[i]));

			float rPhase = (1.0+mu*mu) * (3.0/(16.0*PI));
			const float rhackFactor = 1.0; // <--- XXX HACK XXX
			total += rhackFactor * PI * lightIntensity[i] * lightDiffuse[i] * rInScatter * rPhase;

			// taken from Bruneton (loc. cit.):
			float mPhase = (3.0/(8.0*PI)) * ( (1.0-g*g)*(1.0+mu*mu) ) / ( (2.0+g*g)*pow(1.0+g*g-2.0*g*mu, 1.5));
			const float mhackFactor = 1.0; // <--- XXX HACK XXX
			total += mhackFactor * PI * lightIntensity[i] * lightDiffuse[i] * mInScatter * mPhase;

			total += PI * secondaryScatter * lightIntensity[i] * lightDiffuse[i] * (rInScatter+mInScatter) * constantPhase;
	//}
	

	/*
	// estimate integral of scattering along the line from the eye to the
	// vertex
	mat4 rScatterInt = mat4(0.0);
	mat4 mScatterInt = mat4(0.0);
	vec4 extraIn = vec4(0.0);
	float rScatAtmosInt = 0.0;
	float mScatAtmosInt = 0.0;

	vec4 a,b;
	{
		float atmosNear;
		float atmosFar;
		sphereEntryExitDist(atmosNear, atmosFar, geosphereCenter, eyepos, geosphereAtmosTopRad);
		atmosFar = min(atmosFar, length(eyepos));

		vec3 eyedir = normalize(eyepos);
		vec3 a0 = (atmosFar*eyedir - geosphereCenter) / geosphereRadius;
		vec3 b0 = (atmosNear*eyedir - geosphereCenter) / geosphereRadius;
		a = vec4(a0.x,a0.y,a0.z,0.0);
		b = vec4(b0.x,b0.y,b0.z,0.0);
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
	const int N = 4;
	vec4 d = (b-a)/float(N);
	float len = length(d);

	for (int j=0; j<N+1; j++) {
		float simpson = (j==0||j==N) ? 1.0 : (j == j/2*2) ? 2.0 : 4.0;
		vec4 p = a+float(j)*d;
		float r = length(p);
		if (r < 0.99)
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

		vec4 secondaryLightIntensity = 0.0;
		if (useSecondary == 1)
			secondaryLightIntensity = clamp(d / (2.0*(lightDiscRadii+2.0/rADF)) + 0.5, 0.0, 1.0);

		if (occultedLight == 0) {
			float dist = length(p - occultCentre - y[0]*lightDir[0] );
			lightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
				secondaryLightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
							clamp( ( trad-dist + 2.0/rADF ) / ( trad-absdiff + 4.0/rADF ), 0.0, 1.0)));
		}
		else if (occultedLight == 1) {
			float dist = length(p - occultCentre - y[1]*lightDir[1] );
			lightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
				secondaryLightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
							clamp( ( trad-dist + 2.0/rADF ) / ( trad-absdiff + 4.0/rADF ), 0.0, 1.0)));
		}
		else if (occultedLight == 2) {
			float dist = length(p - occultCentre - y[2]*lightDir[2] );
			lightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
				secondaryLightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
							clamp( ( trad-dist + 2.0/rADF ) / ( trad-absdiff + 4.0/rADF ), 0.0, 1.0)));
		}
		else if (occultedLight == 3) {
			float dist = length(p - occultCentre - y[3]*lightDir[3] );
			lightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
				secondaryLightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
							clamp( ( trad-dist + 2.0/rADF ) / ( trad-absdiff + 4.0/rADF ), 0.0, 1.0)));
		}

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
			// XXX: Entirely ad hoc secondary scatter: split into two components, from below and
			// from above, and in each case just integrate along a line, using density 1/ADF away as
			// estimate of average... there's no good theoretical reasoning behind this, but it
			// seems to work!  FIXME: should include Mie for attenuation at least
			vec4 secondaryScatter = vec4(0.0);
			float samp1 = max(0.0,r-1.0/rADF-1.0);
			float samp2 = r+1.0/rADF;
			float rse1 = exp(-(rADF*samp1));
			float rse2 = exp(-(rADF*samp2));
			float mse1 = exp(-(mADF*samp1));
			float mse2 = exp(-(mADF*samp2));
			secondaryScatter += exp(-(rse1*(r-1)*rextinction + mse1*(r-1)*mextinction)) * (rc * re + mc * me) * (r-1.0);
			float outer = 1.0 + 6.0/rADF;
			if (r < outer)
				secondaryScatter += exp(-(rse2*(outer-r)*rextinction + mse2*(r-1)*mextinction)) * (rc * re + mc * me) * (outer-r);

			extraIn += simpson * (len/3.0) * re * rc * secondaryScatter * (matrixCompMult(attenuation, lightDiffuse) * secondaryLightIntensity);
		}

		if (j==N) {
			vec4 direct = lightIntensity * max(0.0, vec4(tnorm,0.0) * lightDir);
			extraIn += gl_Color * (direct[0] + secondaryLightIntensity[0] * secondaryScatter) * attenuation[0] * lightDiffuse[0];
			extraIn += gl_Color * (direct[1] + secondaryLightIntensity[1] * secondaryScatter) * attenuation[1] * lightDiffuse[1];
			extraIn += gl_Color * (direct[2] + secondaryLightIntensity[2] * secondaryScatter) * attenuation[2] * lightDiffuse[2];
			extraIn += gl_Color * (direct[3] + secondaryLightIntensity[3] * secondaryScatter) * attenuation[3] * lightDiffuse[3];
		}
	}

	vec4 rCol = 0.0;
	mat4 mCol;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		// Actual phase function for Rayleigh scattering is (1 + cos(viewAngle) * (3/(16\pi)), but
		// due to multiple scattering, the dependence on angle is barely present in real skies. As
		// part of our implementation of multiple scattering, therefore, we use the constant phase
		// factor:
		const float rPhase = 1.0/(4.0*PI);
		// We also just multiply up the value by an ad hoc factor. We could try to excuse this by
		// saying it accounts for multiple scatter, but I've no reason other than subjective
		// judgement of results to think it physically realistic!
		const float rhackFactor = 5.0; // <--- XXX HACK XXX
		rCol += rc * rhackFactor * gl_LightSource[i].diffuse * rPhase * rScatterInt[i] * len/3.0;

		// Mie scattering, meanwhile, is highly direction-dependent, so we
		// calculate the phase function in the fragment shader.
		const float mhackFactor = 4.0; // <--- XXX HACK XXX
		mCol[i] = mc * mhackFactor * gl_LightSource[i].diffuse * mScatterInt[i] * len/3.0;
	}

	// For the Earth's atmosphere:
	const float g = 0.76;
	// NOTE: might be worth borrowing from Riley et al, who use the
	// iterated self-convolution of the phase factor to represent multiple
	// forward scattering.

	vec4 atmosDiffuse = rCol;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		vec3 lightDir = normalize(vec3(gl_LightSource[i].position) - geosphereCenter);
		vec3 eyedir = normalize(eyepos);
		float mu = dot(eyedir,lightDir);
		// taken from Bruneton (loc. cit.):
		float mPhase = (3.0/(8.0*PI)) * ( (1-g*g)*(1+mu*mu) ) / ( (2+g*g)*pow(1+g*g-2*g*mu, 1.5));

		atmosDiffuse += mCol[i] * mPhase;
	}
	atmosDiffuse += extraIn;
	const float lightHackFactor = 1.0; // <--- XXX HACK XXX
	atmosDiffuse *= lightHackFactor;
	atmosDiffuse.a = 1.0;
	//float sun = max(0.0, dot(normalize(eyepos),normalize(vec3(gl_LightSource[0].position))));
	gl_FragColor = atmosDiffuse;
	*/

	//vec4 atmosDiffuse = vec4(0.0,0.0,0.0,1.0);
	//{
		//vec3 surfaceNorm = normalize(eyepos - geosphereCenter);
		//for (int i=0; i<NUM_LIGHTS; ++i) {
			//atmosDiffuse += gl_LightSource[i].diffuse * max(0.0, dot(surfaceNorm, normalize(vec3(gl_LightSource[i].position))));
		//}
	//}
	//atmosDiffuse.a = 1.0;
//	float sun = dot(normalize(eyepos),normalize(vec3(gl_LightSource[0].position)));
	gl_FragColor = total + gl_LightModel.ambient*gl_Color +
		//(1.0-fogFactor)*(atmosDiffuse*atmosColor) +
		gl_FrontMaterial.emission;

#ifdef ZHACK
	SetFragDepth(gl_TexCoord[0].z);
#endif
}
