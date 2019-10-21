#include "util.h"


double Util::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle(faceID);
    std::vector <int> pointID;
    VectorT <float,3> points[3];

    for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(fh); curVert.is_valid(); curVert ++)
    {
        pointID.push_back((*curVert).idx());
    }

    for (size_t i=0; i< pointID.size();i++)
    {
        VertexHandle vh = _mesh->vertex_handle(pointID.at(i));
        points[i][0] = _mesh->point(vh)[0];
        points[i][1] = _mesh->point(vh)[1];
        points[i][2] = _mesh->point(vh)[2];
    }
    VectorT<float,3> BmoinsA = points[1] - points[0];
    VectorT<float,3> CmoinsA = points[2] - points[0];

    VectorT<float,3> res;
    res[0] = BmoinsA[1] * CmoinsA[2] - BmoinsA[2] * CmoinsA[1];
    res[1] = BmoinsA[2] * CmoinsA[0] - BmoinsA[0] * CmoinsA[2];
    res[2] = BmoinsA[0] * CmoinsA[1] - BmoinsA[1] * CmoinsA[0];

    return res.norm()/2.0;
}

float Util::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
{
     if(faceID0 == faceID1)
         return 0.0;
     if(vertID0 == vertID1)
         return 0.0;

     VertexHandle vh0 = _mesh->vertex_handle (vertID0);
     FaceHandle fh0 = _mesh->face_handle (faceID0);
     FaceHandle fh1 = _mesh->face_handle (faceID1);
     int Signe = 0;

     MyMesh::FaceVertexCWIter fh = _mesh->fv_cwiter (fh0);

     while ( fh.is_valid() && *fh != vh0 )
         fh++;
     fh++;
     VertexHandle Suivant = *fh;

     if (Suivant.idx () == vertID1)
         Signe = -1;
     else
         Signe = 1;

     OpenMesh::Vec3f normal0 (_mesh->normal(fh0));
     OpenMesh::Vec3f normal1 (_mesh->normal(fh1));

     float scalar = normal0 | normal1;
     return Signe * acos ( scalar );
}

VectorT <float,6> Util::VecteurDirecteursTriangle(MyMesh *_mesh, int vertexID, int faceID){
    FaceHandle fh = _mesh->face_handle(faceID);
    std::vector <int> pointID;
    VectorT <float,3> points[3];
    int A = -1;
    int C = -1, B = -1;

     // On donne au point A son ID
    for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(fh); curVert.is_valid(); curVert ++)
    {
       pointID.push_back((*curVert).idx());
       if(vertexID == (*curVert).idx())
           A = pointID.size()-1;
    }

    // affectation des points trouvés
    for (size_t i=0; i<pointID.size();i++)
    {
        VertexHandle vh = _mesh->vertex_handle(pointID.at(i));
        points[i][0] = _mesh->point(vh)[0];
        points[i][1] = _mesh->point(vh)[1];
        points[i][2] = _mesh->point(vh)[2];
    }

    // On donne aux points B et C leur ID
    for(int i = 0; i<3 ; i++){
        if(i != A){
            if(B == -1)
                B = i;
            else if(C == -1)
                C = i;
        }
    }

    // Calcul des vecteurs AB et AC
    VectorT <float,3> AB = points[B] - points[A];
    VectorT <float,3> AC = points[C] - points[A];
    VectorT <float,6> Vec;
    Vec[0] = AB[0]; Vec[1] = AB[1]; Vec[2] = AB[2];
    Vec[3] = AC[0]; Vec[4] = AC[1]; Vec[5] = AC[2];

    return Vec;
}

float Util::fctK(MyMesh* _mesh, int vertexID){
    float firstPart = 1/AireBarycentrique(_mesh, vertexID);
    float secondPart = 2*M_PI - AngleAbs(_mesh, vertexID);
    return firstPart * secondPart;
}

VectorT <float,3> Util::LongueurArc(MyMesh *_mesh, int vertexID, int vertexID2){

    VectorT <float,3> points[2];

    // affectation des points trouvés
    for (int i=0; i<2;i++)
    {
        VertexHandle vh;

        if(i == 0)
            vh = _mesh->vertex_handle(vertexID);
        else if(i == 1)
            vh = _mesh->vertex_handle(vertexID2);

        points[i][0] = _mesh->point(vh)[0];
        points[i][1] = _mesh->point(vh)[1];
        points[i][2] = _mesh->point(vh)[2];
    }

    // Calcul des vecteurs AB et AC
    VectorT <float,3> AB = points[1] - points[0];

    return AB;
}
float Util::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    VectorT <float,6> vec = VecteurDirecteursTriangle(_mesh, vertexID, faceID);
    VectorT <float,3> AB;
        AB[0] = vec[0];
        AB[1] = vec[1];
        AB[2] = vec[2];
    VectorT <float,3> AC;
        AC[0] = vec[3];
        AC[1] = vec[4];
        AC[2] = vec[5];

    AireBarycentrique(_mesh, vertexID);
    return acos(AB.normalized()|AC.normalized());
}

float Util::fctH(MyMesh* _mesh, int vertexID){

    float firstPart = 1/(4.0*AireBarycentrique(_mesh, vertexID));
    float secondPart = 0;

    VertexHandle v_it = _mesh->vertex_handle(vertexID);
    std::vector<VertexHandle> Vertexs;
    for(MyMesh::VertexVertexIter  vv_it = _mesh->vv_iter(v_it); vv_it; ++vv_it) {
        VertexHandle vh = *vv_it;
        Vertexs.push_back(vh);
    }

   int i = 0;
   for(i =0; i<Vertexs.size()-1;i++) {
       FaceHandle fh;
       FaceHandle fh1;

       bool first = false;
       for (MyMesh::VertexFaceIter curVert = _mesh->vf_iter(v_it); curVert.is_valid(); curVert ++)
       {
           for (MyMesh::VertexFaceIter curVert1 = _mesh->vf_iter(Vertexs.at(i)); curVert1.is_valid(); curVert1 ++)
           {
               if((*curVert).idx() == (*curVert1).idx()){

                   if(!first){
                       fh = curVert;
                       fh1 = curVert;
                       first = true;
                   }
                   else
                       fh1 = curVert;
               }
           }
       }

       int vertexEnFace = Vertexs.at(i).idx();

       OpenMesh::Vec3f VecteurDirecteur = _mesh->point (_mesh->vertex_handle ( vertexEnFace ) ) - _mesh->point ( _mesh->vertex_handle(vertexID));

       OpenMesh::Vec3f normal0 ( _mesh->normal ( fh ) );
       OpenMesh::Vec3f normal1 ( _mesh->normal ( fh1 ) );

       if (((normal0 % normal1) | VecteurDirecteur) < 0)
       {
           FaceHandle tmp = fh;
           fh = fh1;
           fh1 = tmp;
       }
        secondPart += (angleFF(_mesh, fh.idx(), fh1.idx(), vertexID, vertexEnFace) * VecteurDirecteur.norm());
   }
    return firstPart * secondPart;
}

float Util::AireBarycentrique(MyMesh* _mesh, int vertexID){
    float aireTotal;
     VertexHandle v_it = _mesh->vertex_handle(vertexID);

     // parcours des faces autour de vertexID
    for(MyMesh::VertexFaceIter  vf_it = _mesh->vf_iter(v_it); vf_it; ++vf_it) {
        FaceHandle fh = *vf_it;
        aireTotal += faceArea(_mesh, fh.idx());
    }
    return 1/3.0*aireTotal;
}
float Util::AngleAbs(MyMesh* _mesh, int vertexID){
    float angleTotal;
     VertexHandle v_it = _mesh->vertex_handle(vertexID);

    // parcours des faces autour de vertexID
    for(MyMesh::VertexFaceIter  vf_it = _mesh->vf_iter(v_it); vf_it; ++vf_it) {
        FaceHandle fh = *vf_it;
        angleTotal += angleEE(_mesh, v_it.idx(), fh.idx());
    }
    return angleTotal;
}



