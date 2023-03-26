#include "precomp.h"

#include "intersectable.h"

#include "trianglemesh.h"

bool TriangleMesh::TriangleIntersect(const Ray& ray, int tridx, SurfaceInteraction& hit) const {
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
    // calculate b1 parameter and test bounds
    auto b1 = dot(tvec, pvec) * inv_det;
    if (b1 < 0.0f || b1 > 1.0f)
        return false;

    // prepare to test V parameter
    auto qvec = cross(tvec, e1);
    // calculate b2 parameter and test bounds
    auto b2 = dot(ray.D, qvec) * inv_det;
    if (b2 < 0.0f || b1 + b2 > 1.0f)
        return false;
    auto b0 = 1.f - b1 - b2;

    // calculate t, ray intersects triangle
    auto t = dot(e2, qvec) * inv_det;
    if (t <= 0.0f || t >= ray.t)
        return false;

    // compute texture coordinates
    float2 uv[3];
    if (!texcoords.empty()) {
        uv[0] = texcoords[indices[tridx + 0].texcoord_index];
        uv[1] = texcoords[indices[tridx + 1].texcoord_index];
        uv[2] = texcoords[indices[tridx + 2].texcoord_index];
    }
    else {
        uv[0] = float2(0, 0);
        uv[1] = float2(1, 0);
        uv[2] = float2(1, 1);
    }
    auto tc = uv[0] * b0 + uv[1] * b1 + uv[2] * b2;

    // Compute triangle partial derivatives
    float3 dpdu, dpdv;
    // Compute deltas for triangle partial derivatives
    float2 duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
    float3 dp02 = v0 - v2, dp12 = v1 - v2;
    float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
    bool degenerateUV = std::abs(determinant) < 1e-8;
    if (!degenerateUV) {
        float invdet = 1 / determinant;
        dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
        dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }
    if (degenerateUV || sqrLength(cross(dpdu, dpdv)) == 0) {
        // Handle zero determinant for triangle partial derivative matrix
        float3 ng = cross(v2 - v0, v1 - v0);
        if (sqrLength(ng) == 0)
            // The triangle is actually degenerate; the intersection is
            // bogus.
            return false;

        CoordinateSystem(normalize(ng), &dpdu, &dpdv);
    }

    hit = SurfaceInteraction(ray.at(t), tc, -ray.D, dpdu, dpdv, this);
    ray.t = t;

    // compute normal
    if (!normals.empty()) {
        // Initialize Triangle shading geometry

        // Compute shading normal _ns_ for triangle
        auto n0 = normals[indices[tridx + 0].normal_index];
        auto n1 = normals[indices[tridx + 1].normal_index];
        auto n2 = normals[indices[tridx + 2].normal_index];
        float3 ns = n0 * b0 + n1 * b1 + n2 * b2;
        if (sqrLength(ns) > 0.f)
            ns = normalize(ns);
        else
            ns = hit.n;

        // Compute shading tangent _ss_ for triangle
        float3 ss = normalize(hit.dpdu);

        // Compute shading bitangent _ts_ for triangle and adjust _ss_
        float3 ts = cross(ss, ns);
        if (sqrLength(ts) > 0.f) {
            ts = normalize(ts);
            ss = cross(ts, ns);
        }
        else
            CoordinateSystem(ns, &ss, &ts);
        hit.SetShadingGeometry(ss, ts, true);
    }

    return true;
}

shared_ptr<TriangleMesh> TriangleMesh::LoadObj(string inputfile, shared_ptr<Material> mat, const mat4& transform,
    bool ignore_normals, string mtl_search_path, bool verbose) {

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

std::shared_ptr<TriangleMesh> TriangleMesh::CreateBackdrop(const float3& origin, const float3& size, float radius, int steps, std::shared_ptr<Material> material) {
    auto width = size[0];
    auto height = size[1];
    auto depth = size[2];

    // 1. start by generating vertices on both sides of the backdrop
    auto vertices = std::vector<float3>();
    auto normals = std::vector<float3>();
    auto texcoords = std::vector<float2>();

    //    left and right vertices are interleaves in vertices vector
    // 1.a very first point is at (+/- width/2, height, 0)
    vertices.push_back(origin + float3(width / 2, height, 0));
    normals.push_back({ 0, 0, -1 });
    texcoords.push_back({ 0, 0 });
    vertices.push_back(origin + float3(-width / 2, height, 0));
    normals.push_back({ 0, 0, -1 });
    texcoords.push_back({ 0, 1 });

    // Important: We don't want the backwall to shared normals with the bevel steps
    // as it can create weird artifacts for small values of _steps_
    // so we are going to add a small quad just before we reach the bevel
    vertices.push_back(origin + float3(width / 2, radius * 1.1f, 0));
    normals.push_back({ 0, 0, -1 });
    texcoords.push_back({ 1, 0 });
    vertices.push_back(origin + float3(-width / 2, radius * 1.1f, 0));
    normals.push_back({ 0, 0, -1 });
    texcoords.push_back({ 1, 1 });


    // 1.b steps are generated by rotating vector(0, 0, radius) (steps-1) times then moving it to its proper place
    //     for simplification, also include first and last points corresponding to angles 0 and pi/2
    auto stepAngle = PI / (2 * steps);
    for (auto i = 0; i <= steps; i++) {
        // compute rotation of (0, 0, 1) clockwise with angle stepAngle*i
        auto zRot = std::cos(stepAngle * i);
        auto yRot = -std::sin(stepAngle * i); // we are rotating clockwise
        // multiply by radius then move (+/- width/2, radius, -radius)
        vertices.push_back(origin + float3(width / 2, yRot * radius + radius, zRot * radius - radius));
        vertices.push_back(origin + float3(-width / 2, yRot * radius + radius, zRot * radius - radius));

        float3 n = normalize(float3(0, -yRot, -zRot));
        //std::cerr << n << std::endl;
        normals.push_back(n);
        normals.push_back(n);

        texcoords.push_back({ 2.f + i, 0 });
        texcoords.push_back({ 2.f + i, 1 });
    }

    // Important: We don't want the floor to have shared normals with the bevel steps
    // so we are going to add a small quad just after the bevel
    vertices.push_back(origin + float3(width / 2, 0, -radius * 1.1f));
    normals.push_back({ 0, 1, 0 });
    texcoords.push_back({ 3.f + steps, 0 });
    vertices.push_back(origin + float3(-width / 2, 0, -radius * 1.1f));
    normals.push_back({ 0, 1, 0 });
    texcoords.push_back({ 3.f + steps, 1 });

    // 1.c last point is at (+/- width/2, 0, -depth)
    vertices.push_back(origin + float3(width / 2, 0, -depth));
    normals.push_back({ 0, 1, 0 });
    texcoords.push_back({ 4.f + steps, 0 });

    vertices.push_back(origin + float3(-width / 2, 0, -depth));
    normals.push_back({ 0, 1, 0 });
    texcoords.push_back({ 4.f + steps, 1 });

    // 2. generate mesh indices
    auto indices = std::vector<index_type>();
    auto numparts = 4 + steps; // backwall + small quad to avoid artifacts + steps + small quad to avoid artifacts + floor
    for (auto i = 0; i < numparts; i++) {
        // left(i) = i*2
        // right(i) = i*2+1

        // add triangle(left(i), left(i+1), right(i))
        indices.push_back({ i * 2 });
        indices.push_back({ (i + 1) * 2 });
        indices.push_back({ i * 2 + 1 });
        // add triangle(left(i+1), right(i+1), right(i))
        indices.push_back({ (i + 1) * 2 });
        indices.push_back({ (i + 1) * 2 + 1 });
        indices.push_back({ i * 2 + 1 });
    }

    return make_shared<TriangleMesh>(indices, vertices, normals, texcoords, material);
}
