#pragma once

#include "trianglemesh.h"


class Bounds {
public:
    Bounds() {
        bmin3[0] = 1e34f;
        bmin3[1] = 1e34f;
        bmin3[2] = 1e34f;

        bmax3[0] = -1e34f;
        bmax3[1] = -1e34f;
        bmax3[2] = -1e34f;
    }

    bool Intersect(const Ray& ray, float& t) const {
        float tmin = 0.0f;
        float tmax = ray.t;

        for (int a = 0; a < 3; a++) {
            auto t0 = fminf(
                (bmin3[a] - ray.O[a]) / ray.D[a],
                (bmax3[a] - ray.O[a]) / ray.D[a]);
            auto t1 = fmaxf(
                (bmin3[a] - ray.O[a]) / ray.D[a],
                (bmax3[a] - ray.O[a]) / ray.D[a]);
            tmin = fmaxf(t0, tmin);
            tmax = fminf(t1, tmax);
            if ((tmax * 1.00000024f) < tmin)
                return false;
        }
        t = tmin;
        return true;
    }

    void Grow(const Bounds& b) {
        bmin3[0] = fminf(bmin3[0], b.bmin3[0]);
        bmin3[1] = fminf(bmin3[1], b.bmin3[1]);
        bmin3[2] = fminf(bmin3[2], b.bmin3[2]);

        bmax3[0] = fmaxf(bmax3[0], b.bmax3[0]);
        bmax3[1] = fmaxf(bmax3[1], b.bmax3[1]);
        bmax3[2] = fmaxf(bmax3[2], b.bmax3[2]);
    }

    void Grow(const float3& f) {
        bmin3[0] = fminf(bmin3[0], f[0]);
        bmin3[1] = fminf(bmin3[1], f[1]);
        bmin3[2] = fminf(bmin3[2], f[2]);

        bmax3[0] = fmaxf(bmax3[0], f[0]);
        bmax3[1] = fmaxf(bmax3[1], f[1]);
        bmax3[2] = fmaxf(bmax3[2], f[2]);
    }

    float Extend(const int axis) const { return bmax3[axis] - bmin3[axis]; }

    int LongestAxis() const
    {
        int a = 0;
        if (Extend(1) > Extend(0)) a = 1;
        if (Extend(2) > Extend(a)) a = 2;
        return a;
    }

    float3 Centroid() const {
        return {
            (bmin3[0] + bmax3[0]) * 0.5f,
            (bmin3[1] + bmax3[1]) * 0.5f,
            (bmin3[2] + bmax3[2]) * 0.5f
        };
    }

    float3 Offset(const float3& p) const {
        float3 o = p - float3(bmin3[0], bmin3[1], bmin3[2]);
        if (bmax3[0] > bmin3[0]) o.x /= bmax3[0] - bmin3[0];
        if (bmax3[1] > bmin3[1]) o.y /= bmax3[1] - bmin3[1];
        if (bmax3[2] > bmin3[2]) o.z /= bmax3[2] - bmin3[2];
        return o;
    }

    float SurfaceArea() const {
        float3 d(Extend(0), Extend(1), Extend(2));
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }

    static Bounds Union(const Bounds& a, const Bounds& b) {
        Bounds o;
        o.bmin3[0] = fminf(a.bmin3[0], b.bmin3[0]);
        o.bmin3[1] = fminf(a.bmin3[1], b.bmin3[1]);
        o.bmin3[2] = fminf(a.bmin3[2], b.bmin3[2]);

        o.bmax3[0] = fmaxf(a.bmax3[0], b.bmax3[0]);
        o.bmax3[1] = fmaxf(a.bmax3[1], b.bmax3[1]);
        o.bmax3[2] = fmaxf(a.bmax3[2], b.bmax3[2]);
        return o;
    }

    float bmin3[3];
    float bmax3[3];
};

class BVHBuildNode {
public:
    void InitLeaf(int first, int n, const Bounds& b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }

    void InitInterior(shared_ptr<BVHBuildNode> l, shared_ptr<BVHBuildNode> r) {
        left = l;
        right = r;
        bounds = Bounds::Union(l->bounds, r->bounds);
        nPrimitives = 0;
    }

    Bounds bounds;
    shared_ptr<BVHBuildNode> left;
    shared_ptr<BVHBuildNode> right;
    int firstPrimOffset, nPrimitives;
};

struct ALIGN( 32 ) BVHNode {
    Bounds bounds;
    int first; // index of first child for interior nodes, or first primitive for leaf nodes
    int count; // primitives count
};

class Primitive {
public:
    Primitive(const float3& v0, const float3& v1, const float3& v2, int idx) : index(idx) {
        bounds.Grow(v0);
        bounds.Grow(v1);
        bounds.Grow(v2);
        centroid = bounds.Centroid();
    }

public:
    int index; // index of first vertex
    Bounds bounds;
    float3 centroid;
};

struct BucketInfo {
    int count = 0;
    Bounds bounds;
};

class BVHTriMesh : public TriangleMesh {
public:
    BVHTriMesh(shared_ptr<TriangleMesh> trimesh, int maxPrimsInNode = 1) : TriangleMesh(trimesh) {
        // generate primitives vector from mesh indices/vertices
        for (auto i = 0; i < indices.size(); i += 3) {
            primitives.push_back(Primitive(
                vertices[indices[i].vertex_index],
                vertices[indices[i + 1].vertex_index],
                vertices[indices[i + 2].vertex_index], i));
        }

        Timer buildTimer;
        int totalNodes = 0;
        auto root = BuildRecursive(0, (int)primitives.size(), maxPrimsInNode, &totalNodes);
        cerr << "BVH with " << totalNodes << " nodes built in " << (buildTimer.elapsed() * 1000) << " ms\n";

        // Flatten BVH
        buildTimer.reset();

        //nodes = new BVHNode[totalNodes + 1];
        void* alignedNodes = MALLOC64((totalNodes + 1) * sizeof(BVHNode));
        nodes = (BVHNode*)alignedNodes;

        int offset = 2; // place first child at idx 2 so that pair of sibling childs fit in the same 64B cache line
        FlattenBVHTree(root.get(), 0, &offset);
        cerr << "BVH flattened in " << (buildTimer.elapsed() * 1000) << " ms\n";
    }

    ~BVHTriMesh() {
        //delete[] nodes;
        FREE64(nodes);
    }

    virtual bool Intersect(const Ray& ray, Hit& hit) const override {
        float dist;
        if (!nodes[0].bounds.Intersect(ray, dist)) {
            return false;
        }
        return RecursiveHit(nodes[0], ray, hit);
    }

protected:
    shared_ptr<BVHBuildNode> BuildRecursive(int start, int end, int maxPrimsInNode, int* totalNodes);
    void FlattenBVHTree(const BVHBuildNode* node, int offset, int* firstChildOffset);

    bool RecursiveHit(const BVHNode& node, const Ray& ray, Hit& rec) const;

    vector<Primitive> primitives;

    BVHNode* nodes;
};

shared_ptr<BVHBuildNode> BVHTriMesh::BuildRecursive(int start, int end, int maxPrimsInNode, int* totalNodes) {
    auto node = make_shared<BVHBuildNode>();
    (*totalNodes)++;

    // Compute bounds of all primitives in BVH node
    Bounds bounds;
    for (auto i = start; i < end; i++)
        bounds.Grow(primitives[i].bounds);

    int nPrimitives = end - start;
    if (nPrimitives == 1) {
        node->InitLeaf(start, nPrimitives, bounds);
        return node;
    }

    // Compute bound of primitive centroids, choose split dimension _dim_
    Bounds centroidBounds;
    for (auto i = start; i < end; i++)
        centroidBounds.Grow(primitives[i].centroid);
    auto axis = centroidBounds.LongestAxis();

    // if all remaining primitives have the same centroid, create a leaf node
    if (centroidBounds.bmin3[axis] == centroidBounds.bmax3[axis]) {
        node->InitLeaf(start, nPrimitives, bounds);
        return node;
    }

    // Partition primitives using approximate SAH
    auto mid = (start + end) / 2;
    if (nPrimitives <= 2) {
        // if we don't have enough primitives, just split into equal sized subsets
        //TODO we don't really need to call nth_element just for two primitives
        std::nth_element(primitives.begin() + start, primitives.begin() + mid, primitives.begin() + end,
            [axis](const Primitive& a, const Primitive& b) {
                return a.centroid[axis] < b.centroid[axis];
            });
    }
    else {
        // Allocate _BucketInfo_ for SAH partition buckets
        constexpr int nBuckets = 12;
        BucketInfo buckets[nBuckets];

        // Initialize _BucketInfo_ for SAH partition buckets
        for (int i = start; i < end; i++) {
            int b = nBuckets * centroidBounds.Offset(primitives[i].centroid)[axis];
            if (b == nBuckets) b = nBuckets - 1;
            buckets[b].count++;
            buckets[b].bounds.Grow(primitives[i].bounds);
        }

        // Compute costs for splitting after each bucket
        float cost[nBuckets - 1];
        for (int i = 0; i < nBuckets - 1; i++) {
            Bounds b0, b1;
            int count0 = 0, count1 = 0;
            for (int j = 0; j <= i; j++) {
                b0.Grow(buckets[j].bounds);
                count0 += buckets[j].count;
            }
            for (int j = i + 1; j < nBuckets; j++) {
                b1.Grow(buckets[j].bounds);
                count1 += buckets[j].count;
            }
            cost[i] = 1 + (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
        }

        // Find the bucket to split at that minimizes SAH metric
        float minCost = cost[0];
        int minCostSplitBucket = 0;
        for (int i = 1; i < nBuckets - 1; i++) {
            if (cost[i] < minCost) {
                minCost = cost[i];
                minCostSplitBucket = i;
            }
        }

        // Either create leaf or split primitives at selected SAH bucket
        float leafCost = nPrimitives;
        if (nPrimitives > maxPrimsInNode || minCost < leafCost) {
            Primitive* pmid = partition(&primitives[start], &primitives[end - 1] + 1,
                [=](const Primitive& pi) {
                    int b = nBuckets * centroidBounds.Offset(pi.centroid)[axis];
                    if (b == nBuckets) b = nBuckets - 1;
                    return b <= minCostSplitBucket;
                });
            mid = pmid - &primitives[0];
        }
        else {
            // Create leaf node
            node->InitLeaf(start, nPrimitives, bounds);
            return node;
        }
    }

    node->InitInterior(BuildRecursive(start, mid, maxPrimsInNode, totalNodes), BuildRecursive(mid, end, maxPrimsInNode, totalNodes));

    return node;
}

void BVHTriMesh::FlattenBVHTree(const BVHBuildNode* node, int offset, int* firstChildOffset) {
    // store node at offset
    // and its children, if any, at firstChildOffset and firstChildOffset+1
    BVHNode* linearNode = &nodes[offset];
    linearNode->bounds = node->bounds;
    if (node->nPrimitives > 0) {
        // Leaf node
        linearNode->first = node->firstPrimOffset;
        linearNode->count = node->nPrimitives;
    }
    else {
        // Interior node
        linearNode->count = 0;
        linearNode->first = (*firstChildOffset);
        (*firstChildOffset) += 2; // reserve space for both children
        FlattenBVHTree(node->left.get(), linearNode->first, firstChildOffset);
        FlattenBVHTree(node->right.get(), linearNode->first + 1, firstChildOffset);
    }
}

bool BVHTriMesh::RecursiveHit(const BVHNode& node, const Ray& ray, Hit& hit) const {
    auto hit_anything = false;

    // is it a leaf node ?
    if (node.count > 0) {

        for (auto i = 0; i < node.count; i++) {
            auto idx = primitives[node.first + i].index;
            if (TriangleIntersect(ray, idx, hit)) {
                hit_anything = true;
            }
        }

        return hit_anything;
    }

    // it's an interior node
    // check ray/aabb intersection with both children nodes
    BVHNode left = nodes[node.first];
    BVHNode right = nodes[node.first + 1];

    float leftDist;
    bool traverseLeft = left.bounds.Intersect(ray, leftDist);
    float rightDist;
    bool traverseRight = right.bounds.Intersect(ray, rightDist);

    // Check if we should traverse right node first
    bool swap;
    if (traverseLeft && traverseRight)
        swap = rightDist < leftDist;
    else if (traverseLeft || traverseRight)
        swap = !traverseLeft;
    else
        return false; // not sure how often this happens

    if (swap) {
        // start with right child
        BVHNode tmp = left;
        left = right;
        right = tmp;
    }

    // at this point we know for sure we need to traverse first
    if (RecursiveHit(left, ray, hit)) {
        hit_anything = true;
    }
    // only traverse second if (traverseLeft && traverseRight)
    if (traverseLeft && traverseRight && RecursiveHit(right, ray, hit)) {
        hit_anything = true;
    }

    return hit_anything;
}