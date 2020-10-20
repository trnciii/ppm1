#pragma once

#include "vec.h"
#include "constant.h"
#include "data.h"

struct Sphere{
	vec3 p;
	double r;
	double area;
	uint32_t mtlID;

	Sphere(vec3 _p, double _r, uint32_t _m)
	:p(_p), r(_r), mtlID(_m), area(4*kPI*_r*_r){}

	double dist(const Ray &ray) const {
		vec3 OP = p - ray.o;
		double b = dot(ray.d, OP);
		double b2 = b*b;
		double c = dot(OP, OP) - r*r;
		
		if(c < b2){
			double t1 = b - sqrt(b2-c);
			if(0<t1)return t1;

			double t2 = b + sqrt(b2-c);
			if(0<t2)return t2;
		}
		return -1;
	}

	void intersect(Intersection* is, const Ray &ray) const {
		double t = this->dist(ray);
		if((0 <t) && (t< is->dist)){
			is->p = ray.o + t*ray.d;
			is->dist = t;
			is->n = (is->p - p)/r;
			is->mtlID = mtlID;
		}
	}
};