#pragma once

#include <iostream>
#include <iomanip>

#include "vec.h"
#include "data.h"
#include "Scene.h"

void printRule(){std::cout <<"----------  ----------  ----------"<<std::endl;}

void printBr(){putchar('\n');}

void print(const char* s1, const char* s2 = "\n"){
	std::cout <<s1 <<s2;
}

void print(const double& v, const char* s = "\n"){
	double exp = log10(fabs(v));
	if(-3<exp && exp<4) printf("%10.6f%s", v, s);
	else printf("%10.2e%s", v, s);
}

void print(const int& v, const char* s = "\n"){
	double exp = log10(fabs(v));
	if(exp<(4)) printf("%4d%s", v, s);
	else printf("%4e%s", v, s);
}

void print(const uint32_t& v, const char* s = "\n"){
	double exp = log10(fabs(v));
	if(exp<(4)) printf("%4d%s", v, s);
	else printf("%4e%s", v, s);
}

void print(const vec3& a, const char* s = "\n"){
	print(a.x, " ");
	print(a.y, " ");
	print(a.z, s);
}

void print(const Ray& a, const char* s = "\n"){
	print(a.o, "  | ");
	print(a.d, s);
}

void print(const Sphere& sph, const char* str = "\n"){
	print(sph.p, " | ");
	print(sph.r, " | ");
	std::cout <<"mtl: " <<std::setw(2) <<sph.mtlID <<str <<std::flush;
}

void print(const Material& m, const char* str = "\n"){
	std::cout <<std::setw(7);
	if(m.type == Material::Type::EMIT) std::cout <<"EMIT";
	if(m.type == Material::Type::LAMBERT) std::cout <<"LAMBERT";
	
	std::cout <<" | ";
	print(m.color, " | ");
	std::cout <<str <<std::flush;
}

void print(const Camera& a, const char* str="\n"){
	print(a.pos, " | ");
	print(a.basis[2], " | ");
	print(a.flen, str);
}

void print(const Scene& scene){
	puts("[ scene ]");
	printRule();

	puts("camera");
	print(scene.camera);
	printBr();

	puts("materials");
	for(int i=0; i<scene.materials.size(); i++){
		printf("[%2d] ", i);
		print(scene.materials[i]);
	}
	printBr();

	puts("spheres");
	for(int i=0; i<scene.spheres.size(); i++){
		printf("[%2d]", i);
		print(scene.spheres[i]);
	}
	printBr();
	
	printRule();
}