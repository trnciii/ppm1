#pragma once

#include "vec.h"
#include "data.h"

struct Camera{
	vec3 pos;
	vec3 basis[3] = {vec3(1,0,0), vec3(0,0,1), vec3(0,-1,0)};
	double flen = 2;

	Ray ray(double x, double y) const{
		vec3 dir = basis[0]*x + basis[1]*y - basis[2]*flen;
		return Ray(pos, normalize(dir));
	}

	void setDir(vec3 dir, vec3 up){
		basis[0] = normalize(cross(dir, up));
		basis[2] = -normalize(dir);
		basis[1] = cross(basis[2], basis[0]);
	}
};