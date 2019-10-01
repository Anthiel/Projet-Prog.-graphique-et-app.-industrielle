#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>

#define myqDebug() qDebug() << fixed << qSetRealNumberPrecision(4)

/* **** début de la partie à compléter **** */

void MainWindow::nb_elements(MyMesh* _mesh)
{
    qDebug() << _mesh->n_vertices() << "sommets,"<< _mesh->n_faces() << "faces";
}

void MainWindow::verification_voisins(MyMesh* _mesh)
{
    std::vector<int> faces_sans_voisin;
    std::vector<int> points_sans_arrete;
    std::vector<int> arrete_sans_face;

    for(MyMesh::FaceIter curF = _mesh->faces_begin() ; curF != _mesh->faces_end() ; curF ++){
        FaceHandle facf = *curF;
        if(_mesh->valence(facf) == 0){
            faces_sans_voisin.push_back(facf.idx());
         }
    }
    if (faces_sans_voisin.size() == 0){
        qDebug() << "toutes les faces ont au moins un voisin";
    }
    else {
        qDebug() << "Les faces sans faces voisines sont : " << faces_sans_voisin;
    }

    for(MyMesh::VertexIter curV = _mesh->vertices_begin() ; curV != _mesh->vertices_end() ; curV ++){
        int nb_arretes = 0;
        VertexHandle verV = *curV;
        for (MyMesh::VertexEdgeIter verEdge = _mesh->ve_iter(verV); verEdge.is_valid(); verEdge ++){
            nb_arretes += 1;
        }
        if(nb_arretes == 0){
            points_sans_arrete.push_back(verV.idx());
        }
    }
    if (points_sans_arrete.size() == 0){
        qDebug() << "tous les points ont au moins une arrete";
    }
    else {
        qDebug() << "Les points sans arrêtes adjacentes sont : " << points_sans_arrete;
    }

    for(MyMesh::EdgeIter curE = _mesh->edges_begin() ; curE != _mesh->edges_end() ; curE ++){
        EdgeHandle edgE = *curE;
        VertexHandle sommet1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(edgE, 1));
        VertexHandle sommet2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(edgE, 1));

        int face_sommet1 = 0;
        int face_sommet2 = 0;
        for (MyMesh::VertexFaceIter verFace = _mesh->vf_iter(sommet1); verFace.is_valid(); verFace ++){
            face_sommet1 += 1;
        }

        for (MyMesh::VertexFaceIter verFace = _mesh->vf_iter(sommet2); verFace.is_valid(); verFace ++){
            face_sommet2 += 1;
        }

        if (face_sommet1 == 0 && face_sommet2 == 0){
            arrete_sans_face.push_back(edgE.idx());
        }
    }
    if (arrete_sans_face.size() == 0){
        qDebug() << "toutes les arretes appartiennent à au moins une face";
    }
    else {
        qDebug() << "Les arrêtes sans faces adjacentes sont : " << arrete_sans_face;
    }
}

double MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle(faceID);
    std::vector <int> pointID;
    VectorT <float,3> points[3];

    for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(fh); curVert.is_valid(); curVert ++)
    {
        pointID.push_back((*curVert).idx());
    }

    for (int i=0; i<pointID.size();i++)
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

float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
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

VectorT <float,6> MainWindow::VecteurDirecteursTriangle(MyMesh *_mesh, int vertexID, int faceID){
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
    for (int i=0; i<pointID.size();i++)
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

MyMesh::Point MainWindow::barycentreForme(MyMesh* _mesh) {

    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;

    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh = *curVert;
        x += _mesh->point(vh)[0];
        y += _mesh->point(vh)[1];
        z += _mesh->point(vh)[2];
    }

    x /= _mesh->n_vertices();
    y /= _mesh->n_vertices();
    z /= _mesh->n_vertices();

    MyMesh::Point bary(x, y, z);

    MyMesh::VertexHandle newhandle = _mesh->add_vertex(bary);
    _mesh->set_color(newhandle, MyMesh::Color(0,255,0));
    _mesh->data(newhandle).thickness = 2;
    qDebug() << "Nombre de points:" << _mesh->n_vertices();
    myqDebug() << "Barycentre X:" << bary[0];
    myqDebug() << "Barycentre Y:" << bary[1];
    myqDebug() << "Barycentre Z:" << bary[2];
    return bary;
}


void MainWindow::boiteEnglobante(MyMesh* _mesh){
    float minx = _mesh->point(_mesh->vertex_handle(0))[0];
    float miny = _mesh->point(_mesh->vertex_handle(0))[1];
    float minz = _mesh->point(_mesh->vertex_handle(0))[2];

    float maxx = _mesh->point(_mesh->vertex_handle(0))[0];
    float maxy = _mesh->point(_mesh->vertex_handle(0))[1];
    float maxz = _mesh->point(_mesh->vertex_handle(0))[2];

    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh = *curVert;
        if(_mesh->point(vh)[0] < minx) minx = _mesh->point(vh)[0];
        if(_mesh->point(vh)[0] > maxx) maxx = _mesh->point(vh)[0];

        if(_mesh->point(vh)[1] < miny) miny = _mesh->point(vh)[1];
        if(_mesh->point(vh)[1] > maxy) maxy = _mesh->point(vh)[1];

        if(_mesh->point(vh)[2] < minz) minz = _mesh->point(vh)[2];
        if(_mesh->point(vh)[2] > maxz) maxz = _mesh->point(vh)[2];
    }

    float sizex = maxx - minx;
    float sizey = maxy - miny;
    float sizez = maxz - minz;

    float centerx = (maxx + minx) / 2;
    float centery = (maxy + miny) / 2;
    float centerz = (maxz + minz) / 2;

    QVector<float> values = {centerx, centery, centerz,
                            sizex, sizey, sizez};

    qDebug() << values;
}


VectorT <float,3> MainWindow::LongueurArc(MyMesh *_mesh, int vertexID, int vertexID2){

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


float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
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


float MainWindow::fctH(MyMesh* _mesh, int vertexID){

    float firstPart = 1/(4.0*AireBarycentrique(_mesh, vertexID));
    float secondPart = 0;

    VertexHandle v_it = _mesh->vertex_handle(vertexID);
    std::vector<VertexHandle> Vertexs;
    for(MyMesh::VertexVertexIter  vv_it = mesh.vv_iter(v_it); vv_it; ++vv_it) {
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


void MainWindow::H_Curv(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh = *curVert;
        float value = fctH(_mesh, vh.idx());
        _mesh->data(vh).value = value;
    }
}

float MainWindow::AireBarycentrique(MyMesh* _mesh, int vertexID){
    float aireTotal;
     VertexHandle v_it = _mesh->vertex_handle(vertexID);

     // parcours des faces autour de vertexID
    for(MyMesh::VertexFaceIter  vf_it = mesh.vf_iter(v_it); vf_it; ++vf_it) {
        FaceHandle fh = *vf_it;
        aireTotal += faceArea(_mesh, fh.idx());
    }
    return 1/3.0*aireTotal;
}

float MainWindow::AngleAbs(MyMesh* _mesh, int vertexID){
    float angleTotal;
     VertexHandle v_it = _mesh->vertex_handle(vertexID);

    // parcours des faces autour de vertexID
    for(MyMesh::VertexFaceIter  vf_it = mesh.vf_iter(v_it); vf_it; ++vf_it) {
        FaceHandle fh = *vf_it;
        angleTotal += angleEE(_mesh, v_it.idx(), fh.idx());
    }
    return angleTotal;
}

float MainWindow::fctK(MyMesh* _mesh, int vertexID){
    float firstPart = 1/AireBarycentrique(_mesh, vertexID);
    float secondPart = 2*M_PI - AngleAbs(_mesh, vertexID);
    return firstPart * secondPart;
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh = *curVert;
        float value = fctK(_mesh, vh.idx());
        _mesh->data(vh).value = value;
    }
}
/* **** fin de la partie à compléter **** */

void MainWindow::frequence_aires(MyMesh* _mesh)
{
    double maxA = faceArea(&mesh, 0);
    double minA = faceArea(&mesh, 0);
    double diff;
    std::vector <int> compt(10, 0);
    std::vector <double> ratio(10, 0);

    for(MyMesh::FaceIter curF = _mesh->faces_begin() ; curF != _mesh->faces_end() ; curF ++){
        if (faceArea(&mesh, (*curF).idx()) > maxA){
            maxA = faceArea(&mesh, (*curF).idx());
        }
        else if (faceArea(&mesh, (*curF).idx()) < minA){
            minA = faceArea(&mesh, (*curF).idx());
        }
    }

    diff = maxA-minA;

    for(MyMesh::FaceIter curF = _mesh->faces_begin() ; curF != _mesh->faces_end() ; curF ++){
        if (faceArea(&mesh, (*curF).idx()) > minA && faceArea(&mesh, (*curF).idx()) <= minA+((10.0*diff)/100.0)){
            compt[0] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((10.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= minA+((20.0*diff)/100.0)){
            compt[1] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((20.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= minA+((30.0*diff)/100.0)){
            compt[2] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((30.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= minA+((40.0*diff)/100.0)){
            compt[3] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((40.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= minA+((50.0*diff)/100.0)){
            compt[4] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((50.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= minA+((60.0*diff)/100.0)){
            compt[5] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((60.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= minA+((70.0*diff)/100.0)){
            compt[6] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((70.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= minA+((80.0*diff)/100.0)){
            compt[7] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((80.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= minA+((90.0*diff)/100.0)){
            compt[8] ++;
        }
        else if (faceArea(&mesh, (*curF).idx()) > minA+((90.0*diff)/100.0) && faceArea(&mesh, (*curF).idx()) <= maxA){
            compt[9] ++;
        }
    }

    for (unsigned int i = 0; i < compt.size() ; i++){
        ratio[i] = (1.0*compt[i])/(1.0*_mesh->n_faces())*100;
        qDebug() << i*10 << "-" << i*10+10 << " aire : " << ratio[i] << "%";
    }
}

void MainWindow:: frequence_voisinage_sommets(MyMesh* _mesh)
{
    unsigned int maxVal = _mesh->valence(_mesh->vertex_handle(0));
    unsigned int minVal = _mesh->valence(_mesh->vertex_handle(0));

    for(MyMesh::VertexIter curV = _mesh->vertices_begin() ; curV != _mesh->vertices_end() ; curV ++){
        if (_mesh->valence((*curV)) < minVal){
            minVal = _mesh->valence((*curV));
        }
        else if (_mesh->valence((*curV)) > maxVal){
            maxVal = _mesh->valence((*curV));
        }
    }

    std::vector <int> compt(maxVal-minVal, 0);
    std::vector <double> ratio(maxVal-minVal, 0);

    for (MyMesh::VertexIter curV = _mesh->vertices_begin() ; curV != _mesh->vertices_end() ; curV ++){
        compt[_mesh->valence((*curV))-minVal] ++;
    }

    for (unsigned int i = 0; i < compt.size() ; i++){
        ratio[i] = (1.0*compt[i])/(1.0*_mesh->n_vertices())*100;
        qDebug() << "sommets à" << i+minVal << "voisins :" << ratio[i] << "%";
    }
}

/* **** début de la partie boutons et IHM **** */
void MainWindow::on_pushButton_bary_clicked()
{
    barycentreForme(&mesh);
    displayMesh(&mesh, true);
}

/*
    Cette fonction est à utiliser UNIQUEMENT avec le fichier testAngleArea.obj
    Elle est appelée par le bouton "Test angles/aires"

    Elle permet de vérifier les fonctions faceArea, angleFF et angleEE.
    Elle doit afficher :

    Aire de la face 0 : 2
    Aire de la face 1 : 2
    Angle entre les faces 0 et 1 : 1.5708
    Angle entre les faces 1 et 0 : -1.5708
    Angle au sommet 1 sur la face 0 : 0.785398
*/

void MainWindow::on_pushButton_angleArea_clicked()
{
    qDebug() << "Aire de la face 0 :" << faceArea(&mesh, 0);
    qDebug() << "Aire de la face 1 :" << faceArea(&mesh, 1);

    qDebug() << "Angle entre les faces 0 et 1 :" << angleFF(&mesh, 0, 1, 1, 2);
    qDebug() << "Angle entre les faces 1 et 0 :" << angleFF(&mesh, 1, 0, 1, 2);

    qDebug() << "Angle au sommet 1 sur la face 0 :" << angleEE(&mesh, 1, 0);
    qDebug() << "Angle au sommet 3 sur la face 1 :" << angleEE(&mesh, 3, 1);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
    nb_elements(&mesh);
    verification_voisins(&mesh);
    frequence_aires(&mesh);
    frequence_voisinage_sommets(&mesh);
}

/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh *_mesh, bool isTemperatureMap, float mapRange)
{
    //QVector<float> bb_values = boiteEnglobante(&mesh);

    GLuint *triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat *triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat *triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if (isTemperatureMap)
    {
        QVector<float> values;

        if (mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size() * 0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt != fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if (_mesh->data(*fvIt).value > 0)
            {
                triCols[3 * i + 0] = 255;
                triCols[3 * i + 1] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 2] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            else
            {
                triCols[3 * i + 2] = 255;
                triCols[3 * i + 1] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 0] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            if (_mesh->data(*fvIt).value > 0)
            {
                triCols[3 * i + 0] = 255;
                triCols[3 * i + 1] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 2] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            else
            {
                triCols[3 * i + 2] = 255;
                triCols[3 * i + 1] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 0] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            if (_mesh->data(*fvIt).value > 0)
            {
                triCols[3 * i + 0] = 255;
                triCols[3 * i + 1] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 2] = 255 - std::min((_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            else
            {
                triCols[3 * i + 2] = 255;
                triCols[3 * i + 1] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
                triCols[3 * i + 0] = 255 - std::min((-_mesh->data(*fvIt).value / range) * 255.0, 255.0);
            }
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }

        /*for(int bb_index = 0; bb_index < sizeof(bb); bb_index += 3){
            triCols[3 * i + 0] = 0;
            triCols[3 * i + 1] = 0;
            triCols[3 * i + 2] = 255;
            triVerts[3 * i + 0] = bb[bb_index];
            triVerts[3 * i + 1] = bb[bb_index + 1];
            triVerts[3 * i + 2] = bb[bb_index + 2];
            triIndiceArray[i] = i;

            i++;
        }*/
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt != fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3 * i + 0] = _mesh->color(*fIt)[0];
            triCols[3 * i + 1] = _mesh->color(*fIt)[1];
            triCols[3 * i + 2] = _mesh->color(*fIt)[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            triCols[3 * i + 0] = _mesh->color(*fIt)[0];
            triCols[3 * i + 1] = _mesh->color(*fIt)[1];
            triCols[3 * i + 2] = _mesh->color(*fIt)[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
            ++fvIt;
            triCols[3 * i + 0] = _mesh->color(*fIt)[0];
            triCols[3 * i + 1] = _mesh->color(*fIt)[1];
            triCols[3 * i + 2] = _mesh->color(*fIt)[2];
            triVerts[3 * i + 0] = _mesh->point(*fvIt)[0];
            triVerts[3 * i + 1] = _mesh->point(*fvIt)[1];
            triVerts[3 * i + 2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }

        /*for(int bb_index = 0; bb_index < sizeof(bb); bb_index += 3){
            triCols[3 * i + 0] = 0;
            triCols[3 * i + 1] = 0;
            triCols[3 * i + 2] = 255;
            triVerts[3 * i + 0] = bb[bb_index];
            triVerts[3 * i + 1] = bb[bb_index + 1];
            triVerts[3 * i + 2] = bb[bb_index + 2];
            triIndiceArray[i] = i;

            i++;
        }*/
    }

    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint *linesIndiceArray = new GLuint[_mesh->n_edges() * 2 + 12 * 2];
    GLfloat *linesCols = new GLfloat[_mesh->n_edges() * 2 * 3 + 12 * 2 * 3];
    GLfloat *linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3 + 12 * 2 * 3];

    QVector<float> bb_lines = {
        -0.5f, 0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,

        -0.5f, 0.5f, 0.5f,
        -0.5f, -0.5f, 0.5f,

        -0.5f, -0.5f, 0.5f,
        0.5f, -0.5f, 0.5f,

        0.5f, -0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,

        -0.5f, 0.5f, -0.5f,
        0.5f, 0.5f, -0.5f,

        -0.5f, 0.5f, -0.5f,
        -0.5f, -0.5f, -0.5f,

        -0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, -0.5f,

        0.5f, -0.5f, -0.5f,
        0.5f, 0.5f, -0.5f,

        -0.5f, 0.5f, 0.5f,
        -0.5f, 0.5f, -0.5f,

        0.5f, 0.5f, 0.5f,
        0.5f, 0.5f, -0.5f,

        -0.5f, -0.5f, 0.5f,
        -0.5f, -0.5f, -0.5f,

        0.5f, -0.5f, 0.5f,
        0.5f, -0.5f, -0.5f,
    };

    i = 0;
    QHash<float, QList<int>> edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if (t > 0)
        {
            if (!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int>> it(edgesIDbyThickness);
    QList<QPair<float, int>> edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for (int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3 * i + 0] = _mesh->point(vh1)[0];
            linesVerts[3 * i + 1] = _mesh->point(vh1)[1];
            linesVerts[3 * i + 2] = _mesh->point(vh1)[2];
            linesCols[3 * i + 0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3 * i + 1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3 * i + 2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3 * i + 0] = _mesh->point(vh2)[0];
            linesVerts[3 * i + 1] = _mesh->point(vh2)[1];
            linesVerts[3 * i + 2] = _mesh->point(vh2)[2];
            linesCols[3 * i + 0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3 * i + 1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3 * i + 2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    /*for( int bb_index = 0; bb_index < bb_lines.size(); bb_index = bb_index + 3 ) {
        linesVerts[3 * i + 0] = (bb_lines[bb_index] + bb_values[0]) * bb_values[3];
        linesVerts[3 * i + 1] = (bb_lines[bb_index + 1] + bb_values[1]) * bb_values[4];
        linesVerts[3 * i + 2] = (bb_lines[bb_index + 2] + bb_values[2]) * bb_values[5];
        linesCols[3 * i + 0] = 0;
        linesCols[3 * i + 1] = 0;
        linesCols[3 * i + 2] = 255;

        i++;
    }*/

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint *pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat *pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat *pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int>> vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if (t > 0)
        {
            if (!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int>> vitt(vertsIDbyThickness);
    QList<QPair<float, int>> vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for (int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3 * i + 0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3 * i + 1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3 * i + 2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3 * i + 0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3 * i + 1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3 * i + 2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}
