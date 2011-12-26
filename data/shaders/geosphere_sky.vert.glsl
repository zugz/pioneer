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

	// TODO: take into account surface density, and other characteristics of
	// the atmosphere.
	vec4 c = vec4(5.8,13.5,33.1,0)*geosphereScale*geosphereRadius/1000000.0;

	//const float coeffs[5] = {-0.0631,-0.418,8.47,-15.4,10.5};
	//const float sumCoeffs = -0.0631-0.418+8.47-15.4+10.5;
	const float silly = 2.0*exp(-0.0631-0.418+8.47-15.4+10.5)/ADF;
	const int SN = 8;
	vec4 d = (b-a)/float(SN);
	float len = length(d);

	gl_TexCoord[2] = vec4(0.0,0.0,0.0,1.0);
	mat4 lightDir;
	for (int i=0; i<NUM_LIGHTS; ++i)
		lightDir[i] = normalize(gl_LightSource[i].position - (geosphereCenter.x,geosphereCenter.y,geosphereCenter.z,0.0));

	float trad = srad+lrad;
	float absdiff = abs(srad-lrad);

	// estimate integral of scattering along the eyeline
	mat4 scatterInt = mat4(0.0);
	vec4 secondaryScatterInt = 0.0;
	float  scatAtmosInt = 0.0;
	for (int j=0; j<SN+1; j++) {
		vec4 p = a+float(j)*d;
		float r = length(p);
		if (r < 1.0)
			continue;

		float e = exp(-ADF*(r-1.0));

		float lenInvSq = 1.0/dot(p,p);

		vec4 sp = normalize(p);
		float sr = length(p);
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
		// We hand-unroll these loops over lights, to avoid compiler problems:
		primaryScatter[0] = exp(-(lightAtmosInt[0]+scatAtmosInt)*c) * e;
		primaryScatter[1] = exp(-(lightAtmosInt[1]+scatAtmosInt)*c) * e;
		primaryScatter[2] = exp(-(lightAtmosInt[2]+scatAtmosInt)*c) * e;
		primaryScatter[3] = exp(-(lightAtmosInt[3]+scatAtmosInt)*c) * e;
		if (j<SN) scatAtmosInt += e * len / 2.0;

		vec4 d = y*lenInvSq + sqrt((1.0-lenInvSq)*(1.0-(y*y*lenInvSq)));
		vec4 lightIntensity = clamp(d / (2.0*lightDiscRadii) + 0.5, 0.0, 1.0);
		vec4 secondaryLightIntensity = clamp(d / (4.0*lightDiscRadii) + 0.5, 0.0, 1.0);

		if (occultedLight == 0) {
			float dist = length(p - occultCentre - y[0]*lightDir[0] );
			lightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
			secondaryLightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( 2.0 * ( trad-absdiff ) ), 0.0, 1.0)));
		}
		else if (occultedLight == 1) {
			float dist = length(p - occultCentre - y[1]*lightDir[1] );
			lightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
			secondaryLightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( 2.0 * ( trad-absdiff ) ), 0.0, 1.0)));
		}
		else if (occultedLight == 2) {
			float dist = length(p - occultCentre - y[2]*lightDir[2] );
			lightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
			secondaryLightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( 2.0 * ( trad-absdiff ) ), 0.0, 1.0)));
		}
		else if (occultedLight == 3) {
			float dist = length(p - occultCentre - y[3]*lightDir[3] );
			lightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
			secondaryLightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
						clamp( ( trad-dist ) / ( 2.0 * ( trad-absdiff ) ), 0.0, 1.0)));
		}

		if (useSecondary == 1) {
			// Complete hack
			vec4 secondaryScatter = vec4(0.0);
			float se1 = exp(min(0.0,-(ADF*(r-1.0)-1.0)));
			float se2 = exp(-(ADF*(r-1.0)+1.0));
			secondaryScatter += se1*exp(-se1*(1.0/ADF)*c)*c*(r-1.0);
			float outer = 1.0 + 6.0/ADF;
			if (r < outer)
				secondaryScatter += se2*exp(-se2*(1.0/ADF)*c)*c*(outer-r);

			secondaryScatterInt += secondaryScatter * secondaryLightIntensity[0] * primaryScatter[0] * gl_LightSource[0].diffuse;
			secondaryScatterInt += secondaryScatter * secondaryLightIntensity[1] * primaryScatter[1] * gl_LightSource[1].diffuse;
			secondaryScatterInt += secondaryScatter * secondaryLightIntensity[2] * primaryScatter[2] * gl_LightSource[2].diffuse;
			secondaryScatterInt += secondaryScatter * secondaryLightIntensity[3] * primaryScatter[3] * gl_LightSource[3].diffuse;
		}

		primaryScatter[0] *= lightIntensity[0];
		primaryScatter[1] *= lightIntensity[1];
		primaryScatter[2] *= lightIntensity[2];
		primaryScatter[3] *= lightIntensity[3];

		scatterInt += ( primaryScatter ) *
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
	float total = 0.0;
	for (int i=0; i<NUM_LIGHTS; ++i) {
		total += scatterInt[i][0];
		gl_TexCoord[2] += c*scatterInt[i]*len/2.0;
	}
	
	gl_TexCoord[3] = vec4(1.0);
	if (total > 0.0)
		// record ratios
		for (int i=0; i<NUM_LIGHTS; ++i)
			gl_TexCoord[3][i] = scatterInt[i][0]*vec4(1.0/total);

	gl_TexCoord[4] = c * secondaryScatterInt * len/2.0;
	gl_TexCoord[4][3] = 1.0;
}
