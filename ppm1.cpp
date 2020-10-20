// -- todo --
// acceralation structure
// BSDFs
// parallelization
// triangle mesh
// validate disc assumption
// structuration of codes (main)

#include <memory>
#include <algorithm>
#include <filesystem>
#include <cassert>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./stbi/stb_image_write.h"

#include "vec.h"
#include "Scene.h"
#include "Random.h"
#include "print.h"

unsigned char tonemap(double c){
	int c_out = 255*pow(c,(1/2.2)) +0.5;
	if(255 < c_out)c_out = 255;
	if(c_out < 0)c_out = 0;
	return c_out&0xff;
}

int writeImage(vec3* color, int w, int h, const char* name){
	unsigned char *tone = new unsigned char[3*w*h];
	for(int i=0; i<w*h; i++){
		tone[3*i  ] = tonemap(color[i].x);
		tone[3*i+1] = tonemap(color[i].y);
		tone[3*i+2] = tonemap(color[i].z);
	}

	int result = stbi_write_png(name, w, h, 3, tone, 3*w);

	delete[] tone;
	return result;
}

double min(vec3 a){
	if(a.x<a.y && a.x<a.z)return a.z;
	else if(a.y<a.z)return a.y;
	else return a.z;
}
double max(vec3 a){
	if(a.x>a.y && a.x>a.z)return a.x;
	else if(a.y>a.z)return a.y;
	else return a.z;
}

vec3 sampleUniformSphere(double u1, double u2, double* p = nullptr){
	u1 = 2*u1 - 1;
	u2 *= 2*kPI;
	double r = sqrt(1-u1*u1);

	if(p)*p = 0.25/kPI;
	return vec3(r*cos(u2), r*sin(u2), u1);
}

vec3 sampleCosinedHemisphere(double u1, double u2, double* p = nullptr){
	u2 *= 2*kPI;
	double r = sqrt(u1);
	double z = sqrt(1-u1);

	if(p)*p = z/kPI;
	return vec3(r*cos(u2), r*sin(u2), z);
}

vec3 offset(vec3 pos, vec3 dir){return pos + (dir*1e-6);}

void tangentspace(vec3 n, vec3 basis[2]){
	int sg =(n.z < 0) ?-1 :1;
	double a = -1.0/(sg+n.z);
	double b = n.x * n.y * a;
	basis[0] = vec3(
		1.0 + sg * n.x*n.x * a,
		sg*b,
		-sg*n.x
	);
	basis[1] = vec3(
		b,
		sg + n.y*n.y * a,
		-n.y
	);
}

Intersection intersect(const Ray& ray, const Scene& scene){
	Intersection is;
		is.dist = kHUGE;
		is.mtlID = scene.environment;

	for(const Sphere& s : scene.spheres)
		s.intersect(&is, ray);

	if(dot(is.n, ray.d)>0){
		is.n *= -1;
	}

	return is;
}

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
	std::string outDir("result");
	if(!( std::filesystem::exists(outDir) && std::filesystem::is_directory(outDir) )){
		std::cout <<"mkdir " <<outDir <<std::endl;
		printBr();
		assert(std::filesystem::create_directory(outDir));
	}
	
	Scene scene;
	initScene(&scene);

	RNG rand;

	// render parameters
	const int width = 512;
	const int height = 512;
	vec3* result = new vec3[width*height];
	vec3* result_emit = new vec3[width*height];

	int pixelInCenter = (width*(height+1))/2;
	int outInterval = 50;

	int nIteration = 1000;
	int nPhoton = 100;
	int nRay = 4;
	double initialRadius = 3;
	double alpha = 0.8;


	std::vector<hitpoint> hitpoints;
	hitpoints.reserve(width*height*nRay);

	// collect hitpoints
	for(int i=0; i<width*height*nRay; i++){
		int idx = i/nRay;
		int x = idx%width;
		int y = idx/width;

		vec3 loc((double)(2*(x+rand.uniform())-width)/height, (double)-(2*(y+rand.uniform())-height)/height, 0);
		Ray view = scene.camera.ray(loc.x, loc.y);
		Intersection is = intersect(view, scene);
		double pixelFilter = (double)1/nRay;

		if(scene.materials[is.mtlID].type != Material::Type::EMIT)
			hitpoints.push_back(hitpoint(is, initialRadius, pixelFilter, idx, view));
		else result_emit[idx] += scene.materials[is.mtlID].color * pixelFilter;
	}

	// progressive estimation pass
	for(int iteration=1; iteration<=nIteration; iteration++){

		if(iteration%outInterval == 0){
			std::cout <<"itr = " <<iteration <<std::endl;
			printBr();
		}


		std::vector<Photon> photonmap;
		photonmap.reserve(10*nPhoton);

		// construct a photon map
		for(int n=0; n<nPhoton; n++){
			Sphere& source = scene.spheres[scene.lights[0]];
			Ray ray(0,0);
			{
				vec3 N = sampleUniformSphere(rand.uniform(), rand.uniform());
				vec3 P = source.p + (source.r + kTINY)*N;
				vec3 tan[2];
				tangentspace(N, tan);

				vec3 hemi = sampleCosinedHemisphere(rand.uniform(), rand.uniform());

				ray.o = offset(P, N);
				ray.d = (N*hemi.z) + (tan[0]*hemi.x) + (tan[1]*hemi.y);
			}
			
			vec3 ph = scene.materials[source.mtlID].color * source.area * kPI / nPhoton;
			// double initialFlux = max(ph);
			double pTerminate = 1;

			while(rand.uniform() < pTerminate){
			// for(int depth=0; depth<5; depth++){
				Intersection is = intersect(ray, scene);
				Material& mtl = scene.materials[is.mtlID];

				photonmap.push_back(Photon(is.p, ph));

				//update ray and radiant intensity
				{
					vec3 tan[2];
					tangentspace(is.n, tan);
					vec3 hemi = sampleCosinedHemisphere(rand.uniform(), rand.uniform());

					ray.o = offset(is.p, is.n);
					ray.d = (is.n*hemi.z) + (tan[0]*hemi.x) + (tan[1]*hemi.y);
				}
				
				ph *= mtl.color/pTerminate;

				pTerminate = std::min(max(ph), 5.0);
			}
		}

		//accumulate radiance
		for(hitpoint& hit : hitpoints){
			const Material mtl = scene.materials[hit.mtlID];

			int M = 0;
			vec3 tauM = 0;
			for(const Photon& photon : photonmap){
				double d = abs(photon.p, hit.p);
				if(d<hit.R /*&& fabs(dot(hit.n, normalize(photon.p-hit.p))) < hit.R * hit.R * 0.1*/ ){
					double photonFilter = 3*(1 - d/hit.R) / (kPI*hit.R*hit.R); // cone
					// double photonFilter = 1/(kPI*hit.R*hit.R); // constant
					tauM += photonFilter * photon.ph * mtl.color /kPI; // times BSDF
					M++;
				}
			}

			if(hit.N==0){
				hit.N += M;
				hit.tau += tauM;
			}
			else{
				int N = hit.N;
				int Nnext = N + alpha*M;

				double ratio = (double)Nnext/(double)(N+M);

				hit.N = Nnext;
				hit.R *= sqrt(ratio);
				hit.tau = (hit.tau + tauM)*ratio;			
			}

			if( (iteration%outInterval == 0) && (hit.pixel == pixelInCenter) ){
				print(hit.N, " | ");
				print(hit.R, " | ");
				print(hit.tau);
			}
		}


		// output
		if(iteration%outInterval == 0 || iteration == nIteration){

			// std::fill(result, result+(width*height), vec3(0));
			memcpy(result, result_emit, width*height*sizeof(vec3));

			for(hitpoint hit : hitpoints){
				result[hit.pixel] += hit.tau * hit.weight / (iteration);
			}

			printBr();

			char outNum[64];
			sprintf(outNum, "result_%04d.png", iteration);

			std::string out = outDir + "/" + outNum;

			if(writeImage(result, width, height, out.data()) == 1)
				std::cout <<"image saved as " <<out <<std::endl;
			else
				std::cout <<"failed to save image." <<std::endl;


			out = outDir + "/result.png";
			if(writeImage(result, width, height, out.data()) == 1)
				std::cout <<"image saved as " <<out <<std::endl;
			else
				std::cout <<"failed to save image." <<std::endl;
		
			printBr();
			printRule();
		}

	}

	delete[] result;
	delete[] result_emit;
	return 0;
}