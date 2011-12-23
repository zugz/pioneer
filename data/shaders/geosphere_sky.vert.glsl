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
	vec3 a = (skyNear * eyedir - geosphereCenter) / geosphereRadius;
	vec3 b = (skyFar * eyedir - geosphereCenter) / geosphereRadius;

	const float ADF = 500.0;
	float surfaceDensity = atmosColor.w*geosphereAtmosFogDensity;
	vec3 c = vec3(5.8,13.5,33.1)*geosphereScale*geosphereRadius/1000000.0;
	//for (int c=0; c<3; c++)
		//ffac[c] = pow(1.0+float(2-c)/6.0,-4);
	const float silly = 2.0*exp(-0.00287+0.459+3.83-6.80+5.25)/ADF;
	const int SN = 8;
	vec3 d = (b-a)/float(SN);
	float len = length(d);

	gl_TexCoord[2] = vec4(0.0,0.0,0.0,1.0);
	// FIXME: only handling one light!
	for (int i=0; i<1; ++i) {
		vec3 lightDir = normalize(vec3(gl_LightSource[i].position) - geosphereCenter);

		// estimate integral of scattering along the eyeline
		vec3 scatterInt = vec3(0.0,0.0,0.0);
		float scatAtmosInt = 0.0;
		for (int j=0; j<SN+1; j++) {
			vec3 p = a+float(j)*d;
			float r = length(p);
			if (r < 1.0)
				continue;

			float e = exp(-ADF*(r-1.0));
			float y = dot(p,lightDir);
			float a = 1.0 - abs(y)/r;
			float lightAtmosInt = exp(-0.00287 + a*(0.459 + a*(3.83 + a*(-6.80 + a*5.25)))) / ADF;
			if (y < 0.0)
				//lightAtmosInt = 2.0 * exp(-ADF * (sqrt(r*r-y*y)-1.0)) * foo - lightAtmosInt;
				lightAtmosInt = silly - lightAtmosInt;
			lightAtmosInt *= e;

			if (j>0) scatAtmosInt += e * len / 2.0;
			vec3 primaryScatter = exp(-(lightAtmosInt+scatAtmosInt)*c) * e;
			if (j>0) scatAtmosInt += e * len / 2.0;

			float lenInvSq = 1.0/dot(p,p);
			float perp = y;
			float d = perp*lenInvSq + sqrt((1-lenInvSq)*(1-(perp*perp*lenInvSq)));
			float lightIntensity = clamp(d / (2.0*lightDiscRadii[0]) + 0.5, 0.0, 1.0);
			if (occultedLight == i) {
			    vec3 projectedPoint = p - perp*lightDir;
			    float dist = length(projectedPoint - occultCentre);
			    lightIntensity *= (1.0 - mix(0.0, maxOcclusion,
				    clamp( ( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ), 0.0, 1.0)));
			}
			primaryScatter *= lightIntensity;

			vec3 sp = p*(1.0+1.0/ADF);
			float sr = length(sp);
			float se = exp(-ADF*(sr-1.0));
			float sy = dot(sp,lightDir);
			float sa = 1.0 - abs(sy)/sr;
			float secondaryScatterInt = exp(-0.00287 + sa*(0.459 + sa*(3.83 + sa*(-6.80 + sa*5.25)))) / ADF;
			if (sy < 0.0) {
				secondaryScatterInt *= -1;
				secondaryScatterInt += silly;
			}
			secondaryScatterInt *= se;
			vec3 secondaryScatter = exp(-(secondaryScatterInt+scatAtmosInt)*c) * se;
			float slenInvSq = 1.0/dot(sp,sp);
			float sperp = sy;
			float sd = sperp*slenInvSq + sqrt((1-slenInvSq)*(1-(sperp*sperp*slenInvSq)));
			float slightIntensity = clamp(sd / (2.0*lightDiscRadii[0]) + 0.5, 0.0, 1.0);
			if (occultedLight == i) {
			    vec3 projectedPoint = sp - perp*lightDir;
			    float dist = length(projectedPoint - occultCentre);
			    slightIntensity *= (1.0 - mix(0.0, maxOcclusion,
				    clamp( ( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ), 0.0, 1.0)));
			}
			secondaryScatter *= slightIntensity;

			scatterInt += ( primaryScatter +
					secondaryScatter) *
				( (j==0||j==SN) ? 1.0 :
				  2.0 );
		}
		scatterInt *= c*vec3(gl_LightSource[i].diffuse)*len/2.0;
		for (int c=0; c<3; c++)
			gl_TexCoord[2][c] += scatterInt[c];
	}
}
