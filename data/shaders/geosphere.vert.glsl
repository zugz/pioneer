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
}
