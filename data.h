#pragma once

#include "vec.h"

struct Ray;
struct Material;
struct Intersection;
struct hitpoint;
struct Photon;

struct Ray{
	vec3 o;
	vec3 d;

	Ray(vec3 _o, vec3 _d):o(_o), d(_d){}
};

struct Material{

	enum Type{
		EMIT,
		LAMBERT,
		// DIELECTRIC,
		// GGX_REFLECTION,
		// GGX_REFRACTOIN,
	};

	Type type;
	vec3 color = vec3(0.6);
	// double a;
	// double ior;
	
	Material(Material::Type t):type(t){}
};

// ray - object intersection info
struct Intersection{
	float dist;
	vec3 p;
	vec3 n;
	uint32_t mtlID;
	bool backface;
};

// first hit point from eye (following the ppm paper)
struct hitpoint{
	vec3 p;
	vec3 n;
	vec3 wo; // ray direction << done
	uint32_t mtlID;
	uint32_t pixel;
	double R; // current photon radius <<todo
	int N = 0; // accumulated photon count <<todo
	vec3 tau = 0; // accumulated reflected flux <<todo
	double weight;

	hitpoint(Intersection& is, double r, double w, uint32_t px, Ray& ray)
	:p(is.p), n(is.n), wo(-ray.d), mtlID(is.mtlID), pixel(px), N(0), R(r), tau(0), weight(w){}
};

struct Photon{
	vec3 p;
	vec3 ph;

	Photon(vec3 _p, vec3 _ph):p(_p), ph(_ph){}
};