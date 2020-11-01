#pragma once

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#include "vec.h"
#include "Scene.h"
#include "Random.h"
#include "kdtree.h"

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

std::vector<hitpoint> createHitpoints(
	const Scene& scene, const int width, const int height, const int nRay,
	RNG* const rand, vec3* const result_emit)
{
	std::vector<hitpoint> hitpoints;
	hitpoints.reserve(width*height*nRay);
	double initialRadius = 1;

	// collect hitpoints
	for(int i=0; i<width*height*nRay; i++){
		int idx = i/nRay;
		int x = idx%width;
		int y = idx/width;

		vec3 loc((double)(2*(x+rand->uniform())-width)/height, (double)-(2*(y+rand->uniform())-height)/height, 0);
		Ray view = scene.camera.ray(loc.x, loc.y);
		Intersection is = intersect(view, scene);
		double pixelFilter = (double)1/nRay;

		if(scene.materials[is.mtlID].type != Material::Type::EMIT)
			hitpoints.push_back(hitpoint(is, initialRadius, pixelFilter, idx, view));
		else result_emit[idx] += scene.materials[is.mtlID].color * pixelFilter;
	}

	return hitpoints;
}

std::vector<Photon> createPhotonmap(const Scene& scene, int nPhoton, RNG* const rand){
	std::vector<Photon> photonmap;
	photonmap.reserve(10*nPhoton);

	for(int n=0; n<nPhoton; n++){
		const Sphere& source = scene.spheres[scene.lights[0]];
		Ray ray(0,0);
		{
			vec3 N = sampleUniformSphere(rand->uniform(), rand->uniform());
			vec3 P = source.p + (source.r + kTINY)*N;
			vec3 tan[2];
			tangentspace(N, tan);

			vec3 hemi = sampleCosinedHemisphere(rand->uniform(), rand->uniform());

			ray.o = offset(P, N);
			ray.d = (N*hemi.z) + (tan[0]*hemi.x) + (tan[1]*hemi.y);
		}
		
		vec3 ph = scene.materials[source.mtlID].color * source.area * kPI / nPhoton;
		// double initialFlux = max(ph);
		double pTerminate = 1;

		while(rand->uniform() < pTerminate){
		// for(int depth=0; depth<5; depth++){
			Intersection is = intersect(ray, scene);
			const Material& mtl = scene.materials[is.mtlID];

			photonmap.push_back(Photon(is.p, ph, -ray.d));

			//update ray and radiant intensity
			{
				vec3 tan[2];
				tangentspace(is.n, tan);
				vec3 hemi = sampleCosinedHemisphere(rand->uniform(), rand->uniform());

				ray.o = offset(is.p, is.n);
				ray.d = (is.n*hemi.z) + (tan[0]*hemi.x) + (tan[1]*hemi.y);
			}
			
			ph *= mtl.color/pTerminate;

			pTerminate = std::min(max(ph), 5.0);
		}
	}

	return photonmap;
}

void accumulateRadiance(std::vector<hitpoint>& hitpoints, Tree& photonmap, const Scene& scene, const double alpha){
	#pragma omp parallel for schedule(dynamic)
	for(hitpoint& hit : hitpoints){
		const Material& mtl = scene.materials[hit.mtlID];

		int M = 0;
		vec3 tauM = 0;

		std::vector<Photon> candidates = photonmap.searchNN(hit);
		for(const Photon& photon : candidates){
			double d = abs(photon.p, hit.p);
			double photonFilter = 3*(1 - d/hit.R) / (kPI*hit.R*hit.R); // cone
			// double photonFilter = 1/(kPI*hit.R*hit.R); // constant
			tauM += photonFilter * photon.ph * mtl.color /kPI; // times BSDF
			M++;
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
	}
}

// void renderImage