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
    TriangleMesh(
        const vector<index_type>& indices, const vector<float3>& vertices, const vector<float3>& normals,
        const vector<float2>& texcoords, shared_ptr<Material> mat) : indices(move(indices)), vertices(move(vertices)),
        normals(move(normals)), texcoords(move(texcoords)), mat(mat) {};
    TriangleMesh(shared_ptr<TriangleMesh> trimesh) : 
        TriangleMesh(trimesh->indices, trimesh->vertices, trimesh->normals, trimesh->texcoords, trimesh->mat) {};

    virtual bool Intersect(const Ray& ray, Hit& hit) const override {
        bool hit_anything = false;

        for (auto i = 0; i < indices.size(); i+= 3) {
            if (TriangleIntersect(ray, i, hit)) {
                hit_anything = true;
            }
        }

        return hit_anything;
    }

protected:
    bool TriangleIntersect(const Ray& ray, int tridx, Hit& hit) const {
        auto& v0 = vertices[indices[tridx + 0].vertex_index];
        auto& v1 = vertices[indices[tridx + 1].vertex_index];
        auto& v2 = vertices[indices[tridx + 2].vertex_index];

        // find vectors for two edges sharing v0
        auto e1 = v1 - v0;
        auto e2 = v2 - v0;
        // begin calculating determinant - also used to calculate U parameter
        auto pvec = cross(ray.D, e2);
        // if determinant is near zero, ray lies in plane of triangle
        auto det = dot(e1, pvec);

        // check determinant and exit if triangle and ray are parallel
        if (det == 0.0f)
            return false;
        auto inv_det = 1.0f / det;

        // calculate distance from v0 to ray.origin
        auto tvec = ray.O - v0;
        // calculate U parameter and test bounds
        auto u = dot(tvec, pvec) * inv_det;
        if (u < 0.0f || u > 1.0f)
            return false;

        // prepare to test V parameter
        auto qvec = cross(tvec, e1);
        // calculate V parameter and test bounds
        auto v = dot(ray.D, qvec) * inv_det;
        if (v < 0.0f || u + v > 1.0f)
            return false;

        // calculate t, ray intersects triangle
        auto t = dot(e2, qvec) * inv_det;
        if (t <= 0.0f || t >= ray.t)
            return false;

        hit.mat = mat.get();

        // compute texture coordinates
        if (texcoords.empty()) {
            hit.u = u;
            hit.v = v;
        }
        else {
            auto tc0 = texcoords[indices[tridx + 0].texcoord_index];
            auto tc1 = texcoords[indices[tridx + 1].texcoord_index];
            auto tc2 = texcoords[indices[tridx + 2].texcoord_index];
            auto tc = tc0 * (1 - u - v) + tc1 * u + tc2 * v;
            hit.u = tc.x;
            hit.v = tc.y;
        }

        // compute normal
        if (normals.empty()) {
            // mesh doesn't have vertex normals, compute geometric normal
            hit.N = normalize(cross(v1 - v0, v2 - v0));
        }
        else {
            auto n0 = normals[indices[tridx + 0].normal_index];
            auto n1 = normals[indices[tridx + 1].normal_index];
            auto n2 = normals[indices[tridx + 2].normal_index];
            hit.N = normalize(n0 * (1 - u - v) + n1 * u + n2 * v);
        }
        
        hit.I = ray.at(t);
        ray.t = t;

        return true;
    }

    shared_ptr<Material> mat;
    vector<float3> vertices;
    vector<float3> normals;
    vector<float2> texcoords;
    vector<index_type> indices;

public:
    static shared_ptr<TriangleMesh> LoadObj(
        string inputfile,
        shared_ptr<Material> mat,
        const mat4& transform = mat4::Identity(),
        bool ignore_normals = false, 
        string mtl_search_path = "D://models/", 
        bool verbose = false) {

        tinyobj::ObjReaderConfig reader_config;
        reader_config.mtl_search_path = mtl_search_path; // Path to material files

        tinyobj::ObjReader reader;

        if (!reader.ParseFromFile(inputfile, reader_config)) {
            if (!reader.Error().empty()) {
                std::cerr << "TinyObjReader: " << reader.Error();
            }
            exit(1);
        }

        if (!reader.Warning().empty()) {
            cerr << "TinyObjReader: " << reader.Warning();
        }

        auto& attrib = reader.GetAttrib();
        auto& shapes = reader.GetShapes();
        auto& materials = reader.GetMaterials();

        if (verbose) {
            cerr << "attrib.colors.size = " << attrib.colors.size() << endl;
            cerr << "attrib.normals.size = " << attrib.normals.size() << endl;
            cerr << "attrib.texcoords.size = " << attrib.texcoords.size() << endl;
            cerr << "attrib.vertices.size = " << attrib.vertices.size() << endl;
            cerr << "shapes.size = " << shapes.size() << endl;
            cerr << "shapes[0].indices.size = " << shapes[0].mesh.indices.size() << endl;
            cerr << "materials.size = " << materials.size() << endl;
        }

        // extracting vertices and indices
        vector<float3> vertices;
        for (auto i = 0; i < attrib.vertices.size(); i += 3) {
            float3 v(
                attrib.vertices[i + 0], 
                attrib.vertices[i + 1], 
                attrib.vertices[i + 2]
            );
            vertices.push_back(transform.TransformPoint(v));
        }
        vector<index_type> indices;
        for (auto& shape : shapes) {
            for (auto& idx : shape.mesh.indices) {
                indices.push_back({ idx.vertex_index, idx.normal_index, idx.texcoord_index });
            }
        }

        // Normals must be transformed by the inverse transpose of the transformation matrix
        mat4 nTransform = transform.Inverted().Transposed();
        vector<float3> normals;
        if (!ignore_normals) {
            for (auto i = 0; i < attrib.normals.size(); i += 3) {
                float3 n = float3(
                    attrib.normals[i + 0],
                    attrib.normals[i + 1],
                    attrib.normals[i + 2]
                );
                normals.push_back(nTransform.TransformVector(n));
            }
        }

        vector<float2> texcoords;
        for (auto i = 0; i < attrib.texcoords.size(); i += 2) {
            texcoords.push_back(float2(
                attrib.texcoords[i + 0],
                attrib.texcoords[i + 1]
            ));
        }

        return make_shared<TriangleMesh>(indices, vertices, normals, texcoords, mat);
    }
};