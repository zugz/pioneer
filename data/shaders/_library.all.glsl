void DirectionalLight(in int i,
                       in vec3 normal,
                       inout vec4 ambient,
                       inout vec4 diffuse,
                       inout vec4 specular)
{
	float nDotVP; 
	float nDotHV;         
	float pf;             
	nDotVP = max(0.0, dot(normal, normalize(vec3(gl_LightSource[i].position))));
	nDotHV = max(0.0, dot(normal, vec3(gl_LightSource[i].halfVector)));
	if (nDotVP == 0.0) pf = 0.0;
	else pf = pow(nDotHV, gl_FrontMaterial.shininess);
	ambient += gl_LightSource[i].ambient;
	diffuse += gl_LightSource[i].diffuse * nDotVP;
	specular += gl_LightSource[i].specular * pf;
}

void DirectionalLight(in int i,
                       in vec3 normal,
                       inout vec4 ambient,
                       inout vec4 diffuse)
{
	float nDotVP = max(0.0, dot(normal, normalize(vec3(gl_LightSource[i].position))));
	ambient += gl_LightSource[i].ambient;
	diffuse += gl_LightSource[i].diffuse * nDotVP;
}

// Calculate length*density product of a line through the atmosphere
// a - start coord (normalized relative to atmosphere radius)
// b - end coord " "
// centerDensity - atmospheric density at centre of sphere
// length - real length of line in meters
float AtmosLengthDensityProduct(vec3 a, vec3 b, float surfaceDensity, float len)
{
	/* 4 samples */
	float ldprod = 0.0;
	vec3 dir = b-a;
	/* altitude density falloff */
	const float ADF = 500.0;
	// 0.985 = 2.0 - ATMOSPHERE_RADIUS...
	ldprod = surfaceDensity * (
			exp(-ADF*(length(a)-0.985)) +
			exp(-ADF*(length(a + 0.2*dir)-0.985)) +
			exp(-ADF*(length(a + 0.4*dir)-0.985)) +
			exp(-ADF*(length(a + 0.6*dir)-0.985)) +
			exp(-ADF*(length(a + 0.8*dir)-0.985)) +
			exp(-ADF*(length(b)-0.985)));
	ldprod *= len / 6.0;
	return ldprod;	
}

float findSphereEyeRayEntryDistance(in vec3 sphereCenter, in vec3 eyeTo, in float radius)
{
	vec3 v = -sphereCenter;
	vec3 dir = normalize(eyeTo);
	float b = -dot(v, dir);
	float det = (b * b) - dot(v, v) + (radius * radius);
	float entryDist = 0.0;
	if (det > 0.0) {
		det = sqrt(det);
		float i1 = b - det;
		float i2 = b + det;
		if (i2 > 0.0) {
			entryDist = max(i1, 0.0);
		}
	}
	return entryDist;
}


float intensityOfOccultedLight(vec3 lightDir, vec3 v, vec3 occultCentre, float srad, float lrad, float maxOcclusion) {
	vec3 projectedPoint = v - dot(lightDir,v)*lightDir;
	// By our assumptions, the proportion of light blocked at this point by
	// this sphere is the proportion of the disc of radius lrad around
	// projectedPoint covered by the disc of radius srad around occultCentre.
	float dist = length(projectedPoint - occultCentre);

	return 1.0 - mix(0.0, maxOcclusion,
			clamp(
				( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ),
				0.0, 1.0));
}

float intensityOfLightAtPoint(vec3 v, vec3 lightDir, float lightDiscRadius,
		bool eclipse, vec3 occultCentre, float srad, float lrad, float maxOcclusion) {

	// Handle self-shadowing, i.e. "night":
	// d = dot(lightDir,t) where t is the unique point on the unit sphere whose tangent plane
	// contains v, is in the plane of lightDir and d, and is towards the light.
	float lenInvSq = 1.0/dot(v,v);
	float perp = dot(lightDir,v);
	float d = perp*lenInvSq + sqrt((1-lenInvSq)*(1-(perp*perp*lenInvSq)));
	float intensity = 1.0;
	if (lightDiscRadius > 0.0)
		intensity = clamp(d / (2*lightDiscRadius) + 0.5, 0.0, 1.0);

	// Handle eclipse if there is one:
	if (eclipse) {
		// By our assumptions, the proportion of light blocked at this point by
		// this sphere is the proportion of the disc of radius lrad around
		// projectedPoint covered by the disc of radius srad around occultCentre.
		vec3 projectedPoint = v - perp*lightDir;
		float dist = length(projectedPoint - occultCentre);
		intensity *= (1.0 - mix(0.0, maxOcclusion,
					clamp(
						( srad+lrad-dist ) / ( srad+lrad - abs(srad-lrad) ),
						0.0, 1.0)));
	}
	return intensity;
}
