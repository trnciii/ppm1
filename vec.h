#pragma once

#include <math.h>

struct vec3{
	double x;
	double y;
	double z;
	
	vec3():x(0), y(0), z(0) {}
	vec3(double v):x(v), y(v), z(v) {}
	vec3(double _x, double _y, double _z):x(_x), y(_y), z(_z){}

	double operator [] (uint32_t a){return (a==0)? x : (a==1)? y : z;}
};

vec3 operator +(vec3 a, vec3 b){return vec3(a.x+b.x, a.y+b.y, a.z+b.z);}
vec3 operator -(vec3 a, vec3 b){return vec3(a.x-b.x, a.y-b.y, a.z-b.z);}
vec3 operator -(vec3 a){return vec3(-a.x, -a.y, -a.z);}
vec3 operator *(vec3 a, vec3 b){return vec3(a.x*b.x, a.y*b.y, a.z*b.z);}
vec3 operator /(vec3 a, vec3 b){return vec3(a.x/b.x, a.y/b.y, a.z/b.z);}

vec3 operator +=(vec3 &a, vec3 b){return vec3(a.x+=b.x, a.y+=b.y, a.z+=b.z);}
vec3 operator -=(vec3 &a, vec3 b){return vec3(a.x-=b.x, a.y-=b.y, a.z-=b.z);}
vec3 operator *=(vec3 &a, vec3 b){return vec3(a.x*=b.x, a.y*=b.y, a.z*=b.z);}
vec3 operator /=(vec3 &a, vec3 b){return vec3(a.x/=b.x, a.y/=b.y, a.z/=b.z);}

double dot(vec3 a, vec3 b){return a.x*b.x + a.y*b.y + a.z*b.z;}
vec3 cross(vec3 a, vec3 b){return vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);}
double length(vec3 a){return sqrt(a.x*a.x + a.y*a.y +a.z*a.z);}
double length(vec3 a, vec3 b){return length(a-b);}
vec3 normalize(vec3 a){return a/length(a);}