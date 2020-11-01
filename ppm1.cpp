// -- todo --
// acceralation structure
// BSDFs
// parallelization
// triangle mesh
// validate disc assumption
// structuration of codes (main)

#include <filesystem>
#include <cassert>
#include <iostream>

#include "ppm1.h"
#include "print.h"

int initScene(Scene* s){
	s->camera.pos = vec3(0,-10,4);
	s->camera.setDir(vec3(0,1,0), vec3(0,0,1));
	s->camera.flen = 2;

	s->materials[s->environment].type = Material::Type::EMIT;
	s->materials[s->environment].color = vec3(0.05);


	uint32_t light = s->newMaterial(Material::Type::EMIT);
	s->materials[light].color = vec3(30);

	uint32_t white = s->newMaterial(Material::Type::LAMBERT);
	s->materials[white].color = vec3(0.6);

	uint32_t red = s->newMaterial(Material::Type::LAMBERT);
	s->materials[red].color = vec3(0.85, 0.1, 0.1);

	uint32_t green = s->newMaterial(Material::Type::LAMBERT);
	s->materials[green].color = vec3(0.1, 0.85, 0.1);


	// box
	s->add(Sphere(vec3(-1e5, 0, 0), 1e5-4, green)); // left
	s->add(Sphere(vec3( 1e5, 0, 0), 1e5-4, red)); // right
	s->add(Sphere(vec3(0, 0, -1e5), 1e5, white)); // bottom
	s->add(Sphere(vec3(0, 0,  1e5), 1e5-8, white)); // top
	s->add(Sphere(vec3(0,  1e5, 0), 1e5-4, white)); // back
	s->add(Sphere(vec3(1.5, 0.5, 1.5), 1.5, white));
	s->add(Sphere(vec3(0,0,6), 0.5, light)); // light

	print(*s);
	return 0;
}

int main(void){
	// create output directory.
	// using string for dir-name to later create output filename
	// because the format of path::c_str depends on OS.
	std::string outDir("result_master");
	if(!( std::filesystem::exists(outDir) && std::filesystem::is_directory(outDir) )){
		std::cout <<"mkdir " <<outDir <<std::endl;
		printBr();
		assert(std::filesystem::create_directory(outDir));
	}
	
	Scene scene;
	initScene(&scene);

	RNG rand;

	const int width = 512;
	const int height = 512;
	vec3* result = new vec3[width*height];
	vec3* result_emit = new vec3[width*height];

	// render parameters
	int nIteration = 1000;
	int outInterval = 20;

	int nPhoton = 100000;
	int nRay = 4;
	double alpha = 0.7;

	std::vector<hitpoint> hitpoints = createHitpoints(scene, width, height, nRay, &rand, result_emit);

	// progressive estimation pass
	for(int iteration=1; iteration<=nIteration; iteration++){
		std::cout <<"itr = " <<iteration <<std::endl;

		std::vector<Photon> photons = createPhotonmap(scene, nPhoton, &rand);
		std::cout <<"cast photons" <<std::endl;
		Tree photonmap;
		photonmap.copyElements(photons.data(), photons.size());
		photonmap.build();
		std::cout <<"create tree" <<std::endl;
		accumulateRadiance(hitpoints, photonmap, scene, alpha);
		std::cout <<"accumulate radiance" <<std::endl;

		// compose an image
		if(iteration%outInterval == 0 || iteration == nIteration){

			// std::fill(result, result+(width*height), vec3(0));
			memcpy(result, result_emit, width*height*sizeof(vec3));

			for(hitpoint hit : hitpoints){
				result[hit.pixel] += hit.tau * hit.weight / (iteration);
			}

			char outNum[64];
			sprintf(outNum, "result_%04d.png", iteration);
			if(writeImage(result, width, height, (outDir + "/" + outNum).data()) == 1
				&& writeImage(result, width, height, (outDir+"/result.png").data()) == 1)
				std::cout <<"image saved." <<std::endl;
		
			printRule();
		}

	}

	delete[] result;
	delete[] result_emit;
	return 0;
}