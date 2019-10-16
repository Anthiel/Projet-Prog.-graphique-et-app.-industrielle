#ifndef UTIL_H
#define UTIL_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <GL/gl.h>
#include <QVector>
#include <QDebug>

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class Util
{
public:
    static VectorT <float,6> VecteurDirecteursTriangle(MyMesh *_mesh, int vertexID, int faceID);
    static VectorT <float,3> LongueurArc(MyMesh *_mesh, int vertexID, int vertexID2);

    static double faceArea(MyMesh* _mesh, int faceID);
    static float angleFF(MyMesh *_mesh, int faceID0, int faceID1, int vertID0, int vertID1);
    static float angleEE(MyMesh* _mesh, int vertexID, int faceID);
    static void H_Curv(MyMesh* _mesh);
    static void K_Curv(MyMesh* _mesh);
    static int PointEnFace(MyMesh *_mesh, int vertexID, int faceID, int faceID2);
    static float AireBarycentrique(MyMesh* _mesh, int vertexID);
    static float AngleAbs(MyMesh* _mesh, int vertexID);
    static float fctK(MyMesh* _mesh, int vertexID);
    static float fctH(MyMesh* _mesh, int vertexID);
};

#endif // UTIL_H
