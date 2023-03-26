#pragma once

#include "tiny_obj_loader.h"

struct index_type {
    int vertex_index;
    int normal_index;
    int texcoord_index;

    index_type(int idx) :vertex_index(idx), normal_index(idx), texcoord_index(idx) {}
    index_type(int vidx, int nidx, int tidx) : vertex_index(vidx), normal_index(nidx), texcoord_index(tidx) {}
};

class TriangleMesh : public Intersectable {
public:
    TriangleMesh(vector<index_type>& indices, vector<float3>& vertices, vector<float3>& normals, vector<float2>& texcoords, shared_ptr<Material> mat) :
        Intersectable(mat),
        indices(std::move(indices)), 
        vertices(std::move(vertices)), 
        normals(std::move(normals)), 
        texcoords(std::move(texcoords)) {}
    TriangleMesh(shared_ptr<TriangleMesh> trimesh, shared_ptr<Material> mat) : 
        TriangleMesh(trimesh->indices, trimesh->vertices, trimesh->normals, trimesh->texcoords, mat) {};

    virtual bool Intersect(const Ray& ray, SurfaceInteraction& hit) const override {
        bool hit_anything = false;

        for (auto i = 0; i < indices.size(); i+= 3) {
            if (TriangleIntersect(ray, i, hit)) {
                hit_anything = true;
            }
        }

        return hit_anything;
    }

protected:
    bool TriangleIntersect(const Ray& ray, int tridx, SurfaceInteraction& hit) const;

    vector<float3> vertices;
    vector<float3> normals;
    vector<float2> texcoords;
    vector<index_type> indices;

public:
    static shared_ptr<TriangleMesh> LoadObj(string inputfile, shared_ptr<Material> mat, const mat4& transform = mat4::Identity(),
        bool ignore_normals = false, string mtl_search_path = "D://models/", bool verbose = false);

    static std::shared_ptr<TriangleMesh> CreateBackdrop(const float3& origin, const float3& size, float radius, int steps, std::shared_ptr<Material> material);
};