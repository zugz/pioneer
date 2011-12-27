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

	//const float ADF = 500.0;
	float surfaceDensity = atmosColor.w*geosphereAtmosFogDensity;

	// TODO: take into account surface density, and other characteristics of
	// the atmosphere.
	// 
	// Numbers taken from Precomputed Atmospheric Scattering, Bruneton and Neyret 2008
	const float rADF = 795.0;
	const float mADF = 5295.0;
	vec4 rc = vec4(5.8,13.5,33.1,0)*geosphereScale*geosphereRadius/1000000.0;
	vec4 rextinction = rc;
	float mc = 2.0*geosphereScale*geosphereRadius/100000.0;
	float mextinction = mc/0.9;

	//const float coeffs[5] = { 0.10220, -0.81786, 10.53578, -19.373, 12.943};
	//const float mcoeffs[5] = { 1.464499, -3.808514, 25.688546, -48.801948, 29.360771 };
	const vec4 rcoeffs = vec4(-0.81786, 10.53578, -19.373, 12.943);
	const float rconstcoeff = 0.10220;
	const vec4 mcoeffs = vec4( -3.808514, 25.688546, -48.801948, 29.360771 );
	const float mconstcoeff = 1.464499;

	//const float coeffs[5] = {-0.0631,-0.418,8.47,-15.4,10.5};
	//const float sumCoeffs = -0.0631-0.418+8.47-15.4+10.5;
	const float rsilly = 2.0*exp(rconstcoeff+rcoeffs[0]+rcoeffs[1]+rcoeffs[2]+rcoeffs[3])/rADF;
	const float msilly = 2.0*exp(mconstcoeff+mcoeffs[0]+mcoeffs[1]+mcoeffs[2]+mcoeffs[3])/mADF;
	const int SN = 10;
	vec4 d = (b-a)/float(SN);
	float len = length(d);

	mat4 lightDir = mat4(0.0);
	for (int i=0; i<NUM_LIGHTS; ++i)
		lightDir[i] = normalize(gl_LightSource[i].position - (geosphereCenter.x,geosphereCenter.y,geosphereCenter.z,0.0));

	float trad = srad+lrad;
	float absdiff = abs(srad-lrad);

	// estimate integral of scattering along the eyeline
	mat4 rscatterInt = mat4(0.0);
	mat4 mscatterInt = mat4(0.0);
	vec4 secondaryScatterInt = vec4(0.0);
	float rscatAtmosInt = 0.0;
	float mscatAtmosInt = 0.0;
	for (int j=0; j<SN+1; j++) {
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
		vec4 rlightAtmosInt =
			re * ( rsilly*lt + (1.0-2.0*lt) *
			exp(rconstcoeff + a*(rcoeffs[0] + a*(rcoeffs[1] + a*(rcoeffs[2] + a*rcoeffs[3])))) / rADF);
		vec4 mlightAtmosInt =
			me * ( msilly*lt + (1.0-2.0*lt) *
			exp(mconstcoeff + a*(mcoeffs[0] + a*(mcoeffs[1] + a*(mcoeffs[2] + a*mcoeffs[3])))) / mADF);

		if (j>0) rscatAtmosInt += re * len / 2.0;
		if (j>0) mscatAtmosInt += me * len / 2.0;
		mat4 attenuation;
		// We hand-unroll these loops over lights, to avoid compiler problems:
		attenuation[0] = exp(-( (rlightAtmosInt[0]+rscatAtmosInt)*rextinction + (mlightAtmosInt[0]+mscatAtmosInt)*mextinction));
		attenuation[1] = exp(-( (rlightAtmosInt[1]+rscatAtmosInt)*rextinction + (mlightAtmosInt[1]+mscatAtmosInt)*mextinction));
		attenuation[2] = exp(-( (rlightAtmosInt[2]+rscatAtmosInt)*rextinction + (mlightAtmosInt[2]+mscatAtmosInt)*mextinction));
		attenuation[3] = exp(-( (rlightAtmosInt[3]+rscatAtmosInt)*rextinction + (mlightAtmosInt[3]+mscatAtmosInt)*mextinction));
		if (j<SN) rscatAtmosInt += re * len / 2.0;
		if (j<SN) mscatAtmosInt += me * len / 2.0;

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

			secondaryScatterInt += secondaryScatter * secondaryLightIntensity[0] * attenuation[0] * gl_LightSource[0].diffuse;
			secondaryScatterInt += secondaryScatter * secondaryLightIntensity[1] * attenuation[1] * gl_LightSource[1].diffuse;
			secondaryScatterInt += secondaryScatter * secondaryLightIntensity[2] * attenuation[2] * gl_LightSource[2].diffuse;
			secondaryScatterInt += secondaryScatter * secondaryLightIntensity[3] * attenuation[3] * gl_LightSource[3].diffuse;
		}

		attenuation[0] *= lightIntensity[0];
		attenuation[1] *= lightIntensity[1];
		attenuation[2] *= lightIntensity[2];
		attenuation[3] *= lightIntensity[3];

		mat4 rscatter;
		mat4 mscatter;
		rscatter[0] = attenuation[0] * re;
		rscatter[1] = attenuation[1] * re;
		rscatter[2] = attenuation[2] * re;
		rscatter[3] = attenuation[3] * re;
		mscatter[0] = attenuation[0] * me;
		mscatter[1] = attenuation[1] * me;
		mscatter[2] = attenuation[2] * me;
		mscatter[3] = attenuation[3] * me;
		rscatterInt += ( rscatter ) *
			( (j==0||j==SN) ? 1.0 :
			  2.0 );
		mscatterInt += ( mscatter ) *
			( (j==0||j==SN) ? 1.0 :
			  2.0 );

		/*
		vec4 sy = sp * lightDir;
		vec4 sa = vec4(1.0) - abs(sy)/vec4(sr);
		vec4 slt = vec4(lessThan(sy,0.0));
		vec4 secondaryLightInt = silly*slt + (1.0-2.0*slt) *
			exp(-0.0631 + a*(-0.418 + a*(8.47 + a*(-15.4 + a*10.5)))) / ADF;
		secondaryLightInt *= se;

		secondaryScatter[0] = exp(-(secondaryLightInt[0]+scatAtmosInt)*c) * se;
		secondaryScatter[1] = exp(-(secondaryLightInt[1]+scatAtmosInt)*c) * se;
		secondaryScatter[2] = exp(-(secondaryLightInt[2]+scatAtmosInt)*c) * se;
		secondaryScatter[3] = exp(-(secondaryLightInt[3]+scatAtmosInt)*c) * se;

		vec4 sd = sy*slenInvSq + sqrt((1.0-slenInvSq)*(1.0-(sy*sy*slenInvSq)));
		vec4 slightIntensity = clamp(sd / (2.0*lightDiscRadii) + 0.5, 0.0, 1.0);
		if (occultedLight == 0) {
			float dist = length(sp - occultCentre - sy[0]*lightDir[0] );
			slightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
		}
		else if (occultedLight == 1) {
			float dist = length(sp - occultCentre - sy[1]*lightDir[1] );
			slightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
		}
		else if (occultedLight == 2) {
			float dist = length(sp - occultCentre - sy[2]*lightDir[2] );
			slightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
		}
		else if (occultedLight == 3) {
			float dist = length(sp - occultCentre - sy[3]*lightDir[3] );
			slightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff), 0.0, 1.0)));
		}
		secondaryScatter[0] *= slightIntensity[0];
		secondaryScatter[1] *= slightIntensity[1];
		secondaryScatter[2] *= slightIntensity[2];
		secondaryScatter[3] *= slightIntensity[3];

		secondaryScatter *= clamp((r-1)*ADF/3.0, 0.0, 1.0);
		*/

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

	}
	/*
	float total = 0.0;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		total += scatterInt[i][0];
		gl_TexCoord[2] += scatterInt[i]*len/2.0;
	}
	
	gl_TexCoord[3] = vec4(1.0);
	if (total > 0.0)
		// record ratios
		for (int i=0; i<NUM_LIGHTS; ++i)
			gl_TexCoord[3][i] = scatterInt[i][0]*vec4(1.0/total);
	*/

	rCol = 0.0;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		// Actual phase function for Rayleigh scattering is
		// (1 + // cos(viewAngle) * (3/(16\pi)); but due to multiple
		// scattering, this angle dependence is barely present in real skies.
		// As part of our ad-hoc implementation of multiple scattering,
		// therefore, we use a constant phase factor:
		const float rphase = (1.0+0.5)*(3.0/(16*PI));
		rCol += rc * gl_LightSource[i].diffuse * rphase *
			rscatterInt[i] * len/2.0;

		// Mie scattering, meanwhile, is highly direction-dependent, so we
		// calculate the phase function in the fragment shader.
		mCol[i] = mc * gl_LightSource[i].diffuse * mscatterInt[i] * len/2.0;
	}
	secondaryCol = secondaryScatterInt * len/2.0;
}
