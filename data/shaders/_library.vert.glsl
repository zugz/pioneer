
// From http://www.gamedev.net/community/forums/mod/journal/journal.asp?jn=263350&reply_id=3513134
vec4 logarithmicTransform() 
{
	vec4 vertexPosClip = gl_ModelViewProjectionMatrix * gl_Vertex;
	gl_TexCoord[6] = vertexPosClip;
	return vertexPosClip;
}

/*
vec2 findSphereEyeRayEntryExitDistance(in vec3 sphereCenter, in vec3 eyeTo, in float radius)
{
	vec3 v = -sphereCenter;
	vec3 dir = normalize(eyeTo);
	float b = -dot(v, dir);
	float det = (b * b) - dot(v, v) + (radius * radius);
	vec2 dists = vec2(0.0);
	if (det > 0.0) {
		det = sqrt(det);
		float i1 = b - det;
		float i2 = b + det;
		if (i2 > 0) {
			dists.x = max(i1, 0.0);
			dists.y = i2;
		}
	}
	return dists;
}
*/

// opticalDepth: returns estimate of integral of atmosphere density along line
// vertically upwards to infinity, starting at a point p with height h above
// the planet in planet-radii, and with y = cos( angle between integration
// ray and radial ray through p )
//
// Estimation is fairly ad-hoc, but is quite cheap and gets good results -
// according to my tests, ranging over all points and all ADF in the range
// [50,2000], the worst error is about 10%, and typical error is negligible.
float _opticalDepth(float h, float y) {
	// atmosphere density at radius r is exp^(-ADF*(r-1))
	const float ADF=500.0;
	// remaining consts are calculated from ADF (in the future, we may want
	// ADF to vary from planet to planet, in which case these could become
	// uniforms. Note that in reality, it seems that ADF varies even across
	// the earth, due to temperature variation...)
	//
	// goodsin: sine of angle at which we switch from a linear approximation
	// to the integral to an exponential one. Experimentally determined
	// approximate formula for optimal value in terms of ADF:
	// goodsin = 2*(pow(ADF,-0.5))
	const float goodsin = 0.0894;

	// goodcos = sqrt(1-goodsin*goodsin)
	const float goodcos = 0.996;

	// outer: radius at which we consider atmosphere density to be 0
	const float outer = 1.0+7.0/ADF;

	float r = h-1;
	float x = sqrt(h*h-y*y);
    float xsq = x*x;

    float r0 = x/goodcos;
    float h0 = r0+1;
    float y0 = h0*goodsin;
	float outery = sqrt(outer*outer - xsq);


	float int1;
	if (h > outer)
		return 0.0;
	else if (abs(y) < y0) {
		float a0 = (outer/outery-r0/y0)/(outer-r0);
		float d = y0 - abs(y);
		return d*exp(-ADF*sqrt(xsq + (y+d/2.0)*(y+d/2.0))) +
			(exp(-ADF*(h0)) * (r0/y0 - r0*a0 + (r0+(1.0/ADF))*a0))/ADF;
	}
    else {
		float a = (outer/outery-r/y)/(outer-r);
		return (exp(-ADF*h) * (r/y - r*a + (r+(1.0/ADF))*a))/ADF;
	}
}
float opticalDepth(float h, float y) {
	if (y<0.0)
		return 2.0*_opticalDepth(sqrt(h*h-y*y),0) - _opticalDepth(h,-y);
	else
		return _opticalDepth(h,y);
}

float _opticalDepthOReilly(float h, float y) {
	const float ADF=500.0;
	float a = 1.0 - y/(h+1.0);
	// FIXME: wrong constants! They're for ADF=160.
	return exp(-ADF*h) * exp(-0.00287 + a*(0.459 + a*(3.83 + a*(-6.80 + a*5.25)))) / ADF;
}
float _opticalDepthOReilly0(float h) {
	const float ADF=500.0;
	// FIXME: wrong constants! They're for ADF=160.
	//return exp(-ADF*h) * exp(-0.00287 + 0.459 + 3.83 + -6.80 + 5.25) / ADF;
	return exp(-ADF*h) * 0.030854;
}
float opticalDepthOReilly(float h, float y) {
	float d = _opticalDepthOReilly(h,abs(y));
	if (y<0.0)
		// XXX: using h is wrong - should be sqrt( (h+1)*(h+1) - y*y )
		return 2.0*_opticalDepthOReilly0(h) - d;
	else
		return d;
}
