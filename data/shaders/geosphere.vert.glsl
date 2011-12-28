uniform vec3 geosphereCenter;
uniform float geosphereRadius;

uniform int occultedLight;
uniform vec4 occultCentre;
uniform float srad;
uniform float lrad;
uniform float maxOcclusion;
uniform vec4 lightDiscRadii;

void main(void)
{
#ifdef ZHACK
	gl_Position = logarithmicTransform();
#else
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
#endif


	gl_FrontColor = gl_Color;
	gl_TexCoord[1] = gl_ModelViewMatrix * gl_Vertex;
	vec3 tnorm = gl_NormalMatrix * gl_Normal;
	gl_TexCoord[2] = vec4(tnorm.x, tnorm.y, tnorm.z, 0.0);

	// set gl_TexCoord[3][i] to the effective intensity of light i:
	vec4 gc = vec4(geosphereCenter.x, geosphereCenter.y, geosphereCenter.z, 0.0);
	vec4 v = (gl_TexCoord[1] - gc)/geosphereRadius;
	float lenInvSq = 1.0/(length(v)*length(v));
	mat4 lightDir = mat4(0.0);
	for (int i=0; i<NUM_LIGHTS; ++i)
		lightDir[i] = normalize(gl_LightSource[i].position - (geosphereCenter.x,geosphereCenter.y,geosphereCenter.z,0.0));

	vec4 y = v * lightDir;
	// Handle self-shadowing, i.e. "night".
	// d[i] = dot(lightDir[i],t[i]) where t[i] is the unique point on the unit sphere whose tangent
	// plane contains v, is in the plane of lightDir[i] and d, and is towards the light.
	vec4 d = y*lenInvSq + sqrt((1.0-lenInvSq)*(1.0-(y*y*lenInvSq)));
	vec4 lightIntensity = max(vec4(lessThan(lightDiscRadii, 0.0)), clamp(d / (2.0*lightDiscRadii) + 0.5, 0.0, 1.0));

	float trad = srad+lrad;
	float absdiff = abs(srad-lrad);
	if (occultedLight == 0) {
		float dist = length(v - occultCentre - y[0]*lightDir[0] );
		lightIntensity[0] *= (1.0 - mix(0.0, maxOcclusion,
					clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
	}
	if (occultedLight == 1) {
		float dist = length(v - occultCentre - y[1]*lightDir[1] );
		lightIntensity[1] *= (1.0 - mix(0.0, maxOcclusion,
					clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
	}
	if (occultedLight == 2) {
		float dist = length(v - occultCentre - y[2]*lightDir[2] );
		lightIntensity[2] *= (1.0 - mix(0.0, maxOcclusion,
					clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
	}
	if (occultedLight == 3) {
		float dist = length(v - occultCentre - y[3]*lightDir[3] );
		lightIntensity[3] *= (1.0 - mix(0.0, maxOcclusion,
					clamp( ( trad-dist ) / ( trad-absdiff ), 0.0, 1.0)));
	}

	gl_TexCoord[3] = lightIntensity;
}
