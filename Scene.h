#pragma once

#include <memory>
#include <vector>

#include "vec.h"
#include "constant.h"
#include "data.h"
#include "Camera.h"
#include "Sphere.h"


struct Scene{
	Camera camera;
	const uint32_t environment = 0;
	std::vector<Material> materials;
	std::vector<Sphere> spheres;
	std::vector<uint32_t> lights;

	Scene():materials(1, Material(Material::Type::EMIT)){}

	// return nearest hit point
	// Hit intersect(Ray &ray){
	// 	Hit hit;
	// 	hit.dist = kHUGE;
	// 	hit.mtlID = environment;

	// 	for(Sphere s : spheres){
	// 		s.intersect(&hit, ray);
	// 	}

	// 	if(dot(hit.n, ray.d)>0){
	// 		hit.n *= -1;
	// 	}

	// 	return hit;
	// }

	void add(Sphere s){
		if(materials.size() <= s.mtlID)
			throw "invalid material index";

		if(materials[s.mtlID].type == Material::Type::EMIT)
			lights.push_back(spheres.size());

		spheres.push_back(s);
	}

	uint32_t newMaterial(Material::Type t){
		materials.push_back(Material(t));
		return materials.size()-1;
	}
};