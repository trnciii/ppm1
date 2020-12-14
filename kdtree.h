#pragma once

#include <algorithm>
#include <vector>
#include "vec.h"
#include "data.h"

struct Tree{

	struct Node{
		vec3 min;
		vec3 max;

		std::vector<Photon>::iterator begin;
		uint32_t size;
		uint32_t next;

		Node(const std::vector<Photon>::iterator b,
			const std::vector<Photon>::iterator e)
		:min(b[0].p), max(b[0].p), size(e-b), next(0), begin(b)
		{
			for(int i=1; i<size; i++){
				// min = min(min, begin[i]);
				if(begin[i].p.x < min.x) min.x = begin[i].p.x;
				if(begin[i].p.y < min.y) min.y = begin[i].p.y;
				if(begin[i].p.z < min.z) min.z = begin[i].p.z;
				// max = max(max, begin[i]);
				if(max.x < begin[i].p.x) max.x = begin[i].p.x;
				if(max.y < begin[i].p.y) max.y = begin[i].p.y;
				if(max.z < begin[i].p.z) max.z = begin[i].p.z;
			}
		}

		int axis(){
			vec3 dim = max - min;
			return (dim.y<dim.x && dim.z<dim.x)? 0 : ((dim.z<dim.y)? 1 : 2);
		}

		bool intersect(vec3 p, float r){
			return (min.x - r < p.x && p.x < max.x + r
				&& min.y - r < p.y && p.y < max.y + r
				&& min.z - r < p.z && p.z < max.z + r);
		}
	};

	struct Result{
		float distance;
		Photon photon;

		Result(Photon p, float d):distance(d), photon(p){}

		static bool compareDistance(const Result& a, const Result& b){return a.distance < b.distance;}
	};


private:

	std::vector<Photon> verts;
	std::vector<Node> nodes;
	const uint32_t nElements = 1000;

	void split(const std::vector<Photon>::iterator verts_begin,
		const std::vector<Photon>::iterator verts_end,
		const int axis)
	{
		std::sort(verts_begin, verts_end, [axis](Photon a, Photon b){return a.p[axis] < b.p[axis];});	
		std::vector<Photon>::iterator verts_mid =  verts_begin+(verts_end-verts_begin)/2; // split by count

		uint32_t p0 = nodes.size();
		{
			Node node(verts_begin, verts_mid);
			nodes.push_back(node);
			if(nElements < node.size)split(verts_begin, verts_mid, node.axis());
		}
		
		uint32_t p1 = nodes.size();
		{
			Node node(verts_mid, verts_end);
			nodes.push_back(node);
			if(nElements < node.size)split(verts_mid, verts_end, node.axis());
		}

		uint32_t p2 = nodes.size();
	
		nodes[p0].next = p1 - p0;
		nodes[p1].next = p2 - p1;
	};


public:

	bool build(){
		nodes.clear();

		if(verts.size() < 1) return false;
		
		Node node(verts.begin(), verts.end());
		nodes.push_back(node);
		if(nElements < node.size)split(verts.begin(), verts.end(), node.axis());
		nodes[0].next = nodes.size();
		return true;
	}

	std::vector<Result> searchNN(const hitpoint& hit){
		if(!hasTree()) return searchNN_checkAll(hit);

		std::vector<Result> result;

		auto node = nodes.begin();
		while(node < nodes.end()){
			if(node->intersect(hit.p, hit.R)){
				if(node->size <= nElements)
					for(int i=0; i<node->size; i++){
						Photon& photon = node->begin[i];
						vec3 d = photon.p - hit.p;
						double l = abs(d);
						d /= l;

						if(l < hit.R && dot(hit.n, d) < hit.R*0.01)
							result.push_back(Result(node->begin[i], l));
					}

				node++;
			}
			else node += node->next;
		}

		return result;
	}

	std::vector<Result> searchNN_checkAll(const hitpoint& hit){
		std::vector<Result> result;
		
		for(auto v : verts){
			double d = abs(v.p-hit.p);
			if(d < hit.R)
				result.push_back(Result(v, d));
		}

		return result;
	}

	void copyElements(Photon* const elements, uint32_t size){
		std::vector<Photon> v(elements, elements+size);
		verts.swap(v);
		nodes.clear();
	}

	void addElements(Photon* const elements, uint32_t size){
		std::vector<Photon> v(elements, elements+size);
		verts.reserve(verts.size()+size);
		std::copy(v.begin(), v.end(), back_inserter(verts));
		nodes.clear();
	}

	bool hasTree(){return 0 < nodes.size();}

	std::vector<Node> getNodes(){return nodes;}

	std::vector<Photon> getElements(){return verts;}

};