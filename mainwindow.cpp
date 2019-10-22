#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include <iostream>
#include <fstream>

/***************************** Constructor / Destructor *******************************/

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
    this->setWindowTitle("M2GIG - Projet PGAI");

    footerLabel = new QLabel;
    statusBar()->addWidget(footerLabel);

    ui->menuEdition->setEnabled(false);
    ui->menuInfos->setEnabled(false);
    ui->groupBox->setEnabled(false);
    ui->vertexSelect->setEnabled(false);
    ui->edgeSelect->setEnabled(false);
    ui->faceSelect->setEnabled(false);
}
MainWindow::~MainWindow()
{
    delete ui;
}

/***************************** TP1 *******************************/

// Affichage du nombre de sommets et de faces
QVector<int> MainWindow::displayMeshStats(MyMesh* _mesh)
{
    QVector<int> nb_elements = {_mesh->n_vertices(), _mesh->n_faces()};
    return nb_elements;
}

QVector<std::vector<int>> MainWindow::verificationVoisins(MyMesh* _mesh)
{
    QVector<std::vector<int>> list_sans_voisins = {{}, {}, {}};
    std::vector<int> faces_sans_voisin;
    std::vector<int> points_sans_arrete;
    std::vector<int> arrete_sans_face;

    for(MyMesh::FaceIter curF = _mesh->faces_begin() ; curF != _mesh->faces_end() ; curF ++){
        FaceHandle facf = *curF;
        if(_mesh->valence(facf) == 0){
            faces_sans_voisin.push_back(facf.idx());
         }
    }
    if (faces_sans_voisin.size() != 0){
        list_sans_voisins[0] = faces_sans_voisin;
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
    if (points_sans_arrete.size() != 0){
        list_sans_voisins[1] = points_sans_arrete;
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
    if (arrete_sans_face.size() != 0){
        list_sans_voisins[2] = arrete_sans_face;
    }
    return list_sans_voisins;
}

std::vector <double> MainWindow::frequenceAires(MyMesh* _mesh)
{
    double aireTot = 0;
    double maxA = 0;
    double minA = 0;
    double diff;
    std::vector <int> compt(10, 0);
    std::vector <double> ratio(10, 0);

    for(MyMesh::FaceIter curF = _mesh->faces_begin() ; curF != _mesh->faces_end() ; curF ++){
        aireTot += faceArea(&mesh, (*curF).idx());
        if (faceArea(&mesh, (*curF).idx()) > maxA){
            maxA = faceArea(&mesh, (*curF).idx());
        }
        else if (faceArea(&mesh, (*curF).idx()) < minA){
            minA = faceArea(&mesh, (*curF).idx());
        }
    }

    diff = maxA-minA;

    for(MyMesh::FaceIter curF = _mesh->faces_begin() ; curF != _mesh->faces_end() ; curF ++){
        double curr_aire = faceArea(&mesh, (*curF).idx());

        if (curr_aire >= minA && curr_aire <= minA+((10.0*diff)/100.0)){
            compt[0] ++;
        }
        else if (curr_aire > minA+((10.0*diff)/100.0) && curr_aire <= minA+((20.0*diff)/100.0)){
            compt[1] ++;
        }
        else if (curr_aire > minA+((20.0*diff)/100.0) && curr_aire <= minA+((30.0*diff)/100.0)){
            compt[2] ++;
        }
        else if (curr_aire > minA+((30.0*diff)/100.0) && curr_aire <= minA+((40.0*diff)/100.0)){
            compt[3] ++;
        }
        else if (curr_aire > minA+((40.0*diff)/100.0) && curr_aire <= minA+((50.0*diff)/100.0)){
            compt[4] ++;
        }
        else if (curr_aire > minA+((50.0*diff)/100.0) && curr_aire <= minA+((60.0*diff)/100.0)){
            compt[5] ++;
        }
        else if (curr_aire > minA+((60.0*diff)/100.0) && curr_aire <= minA+((70.0*diff)/100.0)){
            compt[6] ++;
        }
        else if (curr_aire > minA+((70.0*diff)/100.0) && curr_aire <= minA+((80.0*diff)/100.0)){
            compt[7] ++;
        }
        else if (curr_aire > minA+((80.0*diff)/100.0) && curr_aire <= minA+((90.0*diff)/100.0)){
            compt[8] ++;
        }
        else if (curr_aire > minA+((90.0*diff)/100.0) && curr_aire <= maxA){
            compt[9] ++;
        }
    }

    std::ofstream outfile (filename + "Aires.csv");

    for (unsigned int i = 0; i < compt.size() ; i++){
        ratio[i] = ((double)compt[i])/((double)_mesh->n_faces())*100.0;
        outfile << i*10 << "-" << i*10+10 << "," << ratio[i] << std::endl;
    }

    outfile.close();

    ratio.push_back(aireTot);

    return ratio;
}

std::vector<double> MainWindow::frequenceVoisinageSommets(MyMesh* _mesh)
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

    std::vector <int> compt(maxVal-minVal+1, 0);
    std::vector <double> ratio(maxVal-minVal+1, 0);

    for (MyMesh::VertexIter curV = _mesh->vertices_begin() ; curV != _mesh->vertices_end() ; curV ++){
        compt[_mesh->valence((*curV))-minVal] ++;
    }

    std::ofstream outfile (filename + "Valence.csv");

    for (unsigned int i = 0; i < compt.size() ; i++){
        ratio[i] = ((double)compt[i])/((double)_mesh->n_vertices())*100.0;
        outfile << " valence de " << minVal + i << " : " << ratio[i] << std::endl;
    }
    ratio.push_back(minVal);
    outfile.close();

    return ratio;
}

// Calcul de la boîte englobante
QVector<float> MainWindow::boiteEnglobante(MyMesh* _mesh)
{
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

    QVector<float> values = {centerx, centery, centerz, sizex, sizey, sizez};
    qDebug() << "Boîte englobante:" << values;


    return values;
}

// Centre gravité / barycentre
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
    qDebug() << "Barycentre X:" << bary[0];
    qDebug() << "Barycentre Y:" << bary[1];
    qDebug() << "Barycentre Z:" << bary[2];
    return bary;
}

void MainWindow::showSelection(MyMesh* _mesh){
    resetAllColorsAndThickness(_mesh);

        if(faceSelection > -1) {
            FaceHandle faceHandle = _mesh->face_handle(faceSelection);
            _mesh->set_color(faceHandle, MyMesh::Color(50,50,255));

            MyMesh::FaceEdgeIter fh_it = mesh.fe_iter(faceHandle);
            for(; fh_it.is_valid(); ++fh_it) {
                EdgeHandle edgeHandle = *fh_it;
                _mesh->set_color(edgeHandle, MyMesh::Color(0, 0, 255));
                _mesh->data(edgeHandle).thickness = 3;
            }

            MyMesh::FaceVertexIter fh_v = mesh.fv_iter(faceHandle);
            for(; fh_v.is_valid(); ++fh_v) {
                VertexHandle vertexHandle = *fh_v;
                _mesh->set_color(vertexHandle, MyMesh::Color(0, 0, 255));
                _mesh->data(vertexHandle).thickness = 7;
            }
        }

        if(edgeSelection > -1) {
            EdgeHandle edgeHandle = _mesh->edge_handle(edgeSelection);
            _mesh->set_color(edgeHandle, MyMesh::Color(0, 255, 0));
            _mesh->data(edgeHandle).thickness = 3;

            HalfedgeHandle heh0 = _mesh->halfedge_handle(edgeHandle, 0); // la première demi-arête
            HalfedgeHandle heh1 = _mesh->halfedge_handle(edgeHandle, 1); // la seconde demi-arête

            VertexHandle v = _mesh->to_vertex_handle(heh0);
            _mesh->set_color(v, MyMesh::Color(0,255,0));
            _mesh->data(v).thickness = 7;
            v = _mesh->to_vertex_handle(heh1);
            _mesh->set_color(v, MyMesh::Color(0,255,0));
            _mesh->data(v).thickness = 7;
        }

        if(vertexSelection > -1) {
            VertexHandle vertexHandle = _mesh->vertex_handle(vertexSelection);
            _mesh->set_color(vertexHandle, MyMesh::Color(255, 0, 0));
            _mesh->data(vertexHandle).thickness = 7;
        }

        // on affiche le nouveau maillage
        displayMesh(_mesh);
}

void MainWindow::showFaceNormal(MyMesh* _mesh){
    FaceHandle currentFace = _mesh->face_handle(faceSelection);
    MyMesh::Normal normal = _mesh->calc_face_normal(currentFace);
    qDebug() << "x: " << normal[0] << " y: " << normal[1] << " z: " << normal[2];
}
void MainWindow::showVertexNormal(MyMesh* _mesh){
    VertexHandle currentVertex = _mesh->vertex_handle(vertexSelection);
    MyMesh::Normal normal = _mesh->calc_vertex_normal(currentVertex);
    qDebug() << "x: " << normal[0] << " y: " << normal[1] << " z: " << normal[2];
}

void MainWindow::deviationNormales(MyMesh* _mesh){
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle currentVertex = *curVert;
        MyMesh::Normal vertexNormal = _mesh->calc_vertex_normal(currentVertex);
        QVector<MyMesh::Normal> faceNormals;
        for(MyMesh::VertexFaceIter curFace = _mesh->vf_iter(currentVertex); curFace.is_valid(); curFace++){
            faceNormals.push_back(_mesh->calc_face_normal(*curFace));
        }

        float maxDeviation = 0.0f;

        for(int i = 0 ; i < faceNormals.size() ; ++i){
            MyMesh::Normal cur = faceNormals[i];
            float scalar = abs(vertexNormal | cur);
            if(acos(scalar) > maxDeviation){
                maxDeviation = acos(scalar);
            }
        }
        _mesh->data(currentVertex).thickness = 5;
        _mesh->set_color(currentVertex, MyMesh::Color(trunc(255*maxDeviation), 0, 0));
    }
    displayMesh(_mesh);
}

std::vector<int> MainWindow::anglesDihedres(MyMesh* _mesh){
    std::vector<int> angleCount(36);
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++){
        EdgeHandle current = *curEdge;
        HalfedgeHandle heh0 = _mesh->halfedge_handle(current, 0);
        HalfedgeHandle heh1 = _mesh->halfedge_handle(current, 1);
        FaceHandle f0 = _mesh->face_handle(heh0);
        FaceHandle f1 = _mesh->face_handle(heh1);
        if(_mesh->is_valid_handle(f0)){
            MyMesh::Normal normal0 = _mesh->calc_face_normal(f0);
            if(_mesh->is_valid_handle(f1)){
                MyMesh::Normal normal1 = _mesh->calc_face_normal(f1);
                float scalar = abs(normal0 | normal1);
                float angle = acos(scalar) * 180 / M_PI;
                int index = abs(floor(angle / 10));
                if (isnan(index) || index < 0 ) index = 0;
                angleCount[index] += 1;
            }
            else{
                continue;
            }
        }
        else{
            continue;
        }
    }

    for(unsigned int i = 0 ; i < angleCount.size() ; ++i){
        qDebug() << i << ": " << angleCount[i];
    }

    std::ofstream outfile (filename + "Angles.csv");

    for(unsigned int i = 0 ; i < angleCount.size() ; ++i){
         outfile << i * 10 << "-" << (i + 1) * 10 << "," << angleCount[i] <<  std::endl;
    }

    outfile.close();

    return angleCount;
}
// Calcul moyenne des normales aux faces concourantes


float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle ( faceID );

    std::vector<VertexHandle> vertexes;

    // On récupère les 3 sommets de la face
    MyMesh::FaceVertexIter fh_v = _mesh->fv_iter(fh);
    for(; fh_v.is_valid(); ++fh_v)
        vertexes.push_back ( *fh_v );

    // On créé deux vecteurs avec le même point de départ
    OpenMesh::Vec3f vectorAB = _mesh->point(vertexes[1]) - _mesh->point(vertexes[0]);
    OpenMesh::Vec3f vectorAC = _mesh->point(vertexes[2]) - _mesh->point(vertexes[0]);

    // On calcule le produit vectoriel de ses deux vecteurs et retourne sa norme divisée par 2
    OpenMesh::Vec3f product = vectorAB % vectorAC;
    float norm = product.norm();

    return norm / 2.0f;
}

// Calcule de l'aire barycentrique sous un sommet
float MainWindow::baryArea(MyMesh* _mesh, int vertID){
    float baryArea = 0;

    VertexHandle vh = _mesh->vertex_handle ( vertID );

    // On somme toutes les aires des faces voisines au sommet
    MyMesh::VertexFaceIter vf = _mesh->vf_iter ( vh );
    for ( ; vf.is_valid ( ) ; ++vf ) {
        FaceHandle current = *vf;
        baryArea += faceArea ( _mesh , current.idx( ) );
    }

    // On retourne cette somme divisée par 3
    return baryArea / 3.0f;
}

// Calcul de l'angle entre deux faces
float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
{
    int sign = 0;

    VertexHandle vh0 = _mesh->vertex_handle ( vertID0 );
    FaceHandle fh0 = _mesh->face_handle ( faceID0 );


    // On créé un itérateur clockwise sur les sommets de la face d'handle faceID0
    // et on tourne jusqu'à arriver sur le sommet d'handle vertID0
    MyMesh::FaceVertexCWIter fh_cwv = _mesh->fv_cwiter ( fh0 );
    while ( fh_cwv.is_valid ( ) && *fh_cwv != vh0 ) ++fh_cwv;

    VertexHandle next = *++fh_cwv;

    // Si le suivant du sommet d'handle vertID0 sur la face d'id faceID0 est le sommet
    // d'handle vertID1, on est dans le sens négatif des faces, le signe de l'angle
    // est négatif. Sinon, il est positif.
    if ( next.idx ( ) == vertID1 ) sign = -1;
    else sign = 1;

    // On récupère les normales, calcule leur produit scalaire et on renvoie l'angle
    // signé avec le produit scalaire et le signe trouvé plus tot
    OpenMesh::Vec3f normal0 (_mesh->normal ( fh0 ) );
    OpenMesh::Vec3f normal1 (_mesh->normal ( _mesh->face_handle ( faceID1 ) ) );
    float scalar = normal0 | normal1;

    return sign * acos ( scalar );
}

float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    FaceHandle fh = _mesh->face_handle ( faceID );
    VertexHandle vh = _mesh->vertex_handle ( vertexID );
    std::vector<VertexHandle> vertexes;

    // On récupère les sommets de la faces qui ne sont pas le sommet d'ID vertexID
    MyMesh::FaceVertexIter fh_v = _mesh->fv_iter(fh);
    for(; fh_v.is_valid(); ++fh_v) {
        VertexHandle current = *fh_v;
        if( current.idx() != vertexID )
            vertexes.push_back ( current );
    }

    // On retourne l'angle calculé depuis le produit scalaire
    OpenMesh::Vec3f vectorAB = _mesh->point(vertexes[0]) - _mesh->point(vh);
    OpenMesh::Vec3f vectorAC = _mesh->point(vertexes[1]) - _mesh->point(vh);
    return acos ( vectorAB.normalize() | vectorAC.normalize() );
}

void MainWindow::H_Curv(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++) {
        MyMesh::VertexHandle current = *curVert;
        float val = 0.0f;

        for (MyMesh::VertexEdgeIter currentEdge = _mesh->ve_iter ( current ); currentEdge.is_valid(); currentEdge++)
        {
            MyMesh::EdgeHandle eh = *currentEdge;
            MyMesh::HalfedgeHandle heh0 = _mesh->halfedge_handle(eh, 0);
            MyMesh::HalfedgeHandle heh1 = _mesh->halfedge_handle(eh, 1);

            FaceHandle fh0 = _mesh->face_handle(heh0);
            FaceHandle fh1 = _mesh->face_handle(heh1);

            // Si l'arête est en bordure, on ne traite qu'une face
            if ( fh1.idx ( ) > _mesh->n_faces ( ) )
                fh1 = fh0;

            // On détermine le sommet opposé au sommet courant sur l'arête courante
            int vertex2ID = _mesh->to_vertex_handle(heh1).idx();
            if (vertex2ID == current.idx ( ) )
                vertex2ID = _mesh->to_vertex_handle(heh0).idx();

            // Vecteur entre le sommet courant et son opposé sur l'arête courante
            OpenMesh::Vec3f currentOppVector = _mesh->point ( _mesh->vertex_handle ( vertex2ID ) ) - _mesh->point ( current );

            OpenMesh::Vec3f normal0 ( _mesh->normal ( fh0 ) );
            OpenMesh::Vec3f normal1 ( _mesh->normal ( fh1 ) );

            // On ordonne les faces
            if ( ( ( normal0 % normal1 ) | currentOppVector ) < 0 )
            {
                FaceHandle tempF = fh0;
                fh0 = fh1;
                fh1 = tempF;
            }

            // On somme les valeurs des angles entre les faces * le vecteur entre le sommet et son opposé
            val += currentOppVector.norm ( ) * angleFF ( _mesh , fh0.idx ( ) , fh1.idx ( ) , current.idx() , vertex2ID );
        }

        // On calcule le résultat sur le sommet et on le lui attribue dans les propriétés du maillage
        val /= ( 4 * baryArea ( _mesh , current.idx ( ) ) );
        _mesh->data ( current ).value = val;
    }
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    for ( MyMesh::VertexIter curVert = _mesh->vertices_begin() ; curVert!=_mesh->vertices_end() ; ++curVert ) {
        VertexHandle current = *curVert;

        // On calcule l'angle barycentrique sous le sommet courant
        float area = baryArea ( _mesh , current.idx() );
        float angleEESum = 0;

        //On itère les faces autour du sommet courant et somme les angles entre leurs arêtes
        MyMesh::VertexFaceIter vf = _mesh->vf_iter ( current );
        for ( ; vf.is_valid ( ) ; ++vf ) {
            FaceHandle currentFace = *vf;
            angleEESum += angleEE ( _mesh , current.idx ( ) , currentFace.idx ( ) );
        }

        // On calcule le résultat sur le sommet et on le lui attribue dans les propriétés du maillage
        _mesh->data ( current ).value = ( 1 / area ) * ( 2 * M_PI - angleEESum );
    }
}

/***************************** UTILS *******************************/

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
            std::sort(values.begin(), values.end());
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

/***************************** IHM *******************************/

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));
    if(fileName.isEmpty()) return;
    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());
    filename = fileName.toUtf8().constData();

    mesh.request_vertex_normals();
    mesh.request_face_normals();
    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    this->resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    this->displayMesh(&mesh);

    this->displayMeshStats(&mesh);
    this->verificationVoisins(&mesh);

    QString nb_elements = QString::number(displayMeshStats(&mesh)[0]) + " sommets, " + QString::number(displayMeshStats(&mesh)[1]) + " faces";
    footerLabel->setText(nb_elements);

    ui->menuEdition->setEnabled(true);
    ui->menuInfos->setEnabled(true);
    ui->groupBox->setEnabled(true);
    ui->vertexSelect->setEnabled(true);
    ui->edgeSelect->setEnabled(true);
    ui->faceSelect->setEnabled(true);
}

void MainWindow::on_pushButton_box_clicked()
{
    QVector<float> result = boiteEnglobante(&mesh);
    QString text("Boîte englobante: ");
    for(int i = 0; i < result.size(); i++){
        text.append(QString::number(result[i]) + " ");
    }
    QMessageBox::information(this, tr("Boîte englobante"), text);
}

void MainWindow::on_pushButton_bary_clicked()
{
    MyMesh::Point point = barycentreForme(&mesh);
    QMessageBox::information(this,
                             tr("Barycentre"),
                             QString("Le barycentre est de coordonnées (" + QString::number(point[0]) + ", " + QString::number(point[1]) + ", " + QString::number(point[2]) + ")"));
    this->displayMesh(&mesh, true);
}

void MainWindow::on_pushButton_aire_clicked()
{
    std::vector<double> ratio = frequenceAires(&mesh);
    QString text("Fréquences des aires : \n");
    for(unsigned int i = 0; i < ratio.size()-1; i++){
        text.append("Aire entre " + QString::number(i*10) + "% et " + QString::number(i*10+10) + "% : " + QString::number(ratio[i]) + "% \n");
    }
    QMessageBox::information(this, tr("Fréquence aires"), text);

    ui->plainTextEdit->setPlainText("Aire totale du maillage = " + QString::number(ratio[ratio.size()-1]));
}

void MainWindow::on_pushButton_freq_valence_clicked()
{
    std::vector<double> ratio = frequenceVoisinageSommets(&mesh);
    QString text("Valences des sommets : \n");
    for(unsigned int i = 0; i < ratio.size()-1; i++){
        text.append("Valence de " + QString::number(ratio[ratio.size()-1] +i) + " : " + QString::number(ratio[i]) + "% \n");
    }
    QMessageBox::information(this, tr("Fréquence valences sommets"), text);
}

void MainWindow::on_vertexSelect_valueChanged(int arg1)
{
    vertexSelection = arg1;
    showSelection(&mesh);
}

void MainWindow::on_edgeSelect_valueChanged(int arg1)
{
    edgeSelection = arg1;
    showSelection(&mesh);
    showVertexNormal(&mesh);
}

void MainWindow::on_faceSelect_valueChanged(int arg1)
{
    faceSelection = arg1;
    showSelection(&mesh);
    showFaceNormal(&mesh);
}

void MainWindow::on_pushButton_dev_clicked()
{
     deviationNormales(&mesh);
}

void MainWindow::on_pushButton_dihedral_clicked()
{
    std::vector<int> list_angle = anglesDihedres(&mesh);
    QString text("Angles dièdres : \n");
    for(int i = 0; i < list_angle.size(); i++){
        text.append(QString::number(i) + " : " + QString::number(list_angle[i]) + "\n");
    }
    QMessageBox::information(this, tr("Angles dièdres"), text);
}

void MainWindow::on_pushButton_h_clicked()
{
    H_Curv(&mesh);
    displayMesh(&mesh, true);
}

void MainWindow::on_pushButton_k_clicked()
{
    K_Curv(&mesh);
    displayMesh(&mesh, true);
}


void MainWindow::on_actionCharger_OBJ_triggered()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());
    filename = fileName.toUtf8().constData();

    mesh.request_vertex_normals();
    mesh.request_face_normals();
    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    this->resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    this->displayMesh(&mesh);

    this->displayMeshStats(&mesh);

    QString nb_elements = QString::number(displayMeshStats(&mesh)[0]) + " sommets, " + QString::number(displayMeshStats(&mesh)[1]) + " faces";
    footerLabel->setText(nb_elements);

    ui->menuEdition->setEnabled(true);
    ui->menuInfos->setEnabled(true);
    ui->groupBox->setEnabled(true);
    ui->vertexSelect->setEnabled(true);
    ui->edgeSelect->setEnabled(true);
    ui->faceSelect->setEnabled(true);
}

void MainWindow::on_actionCourbure_Moyenne_triggered()
{
    H_Curv(&mesh);
    displayMesh(&mesh, true);
}

void MainWindow::on_actionCourbure_Gaussienne_triggered()
{
    K_Curv(&mesh);
    displayMesh(&mesh, true);
}

void MainWindow::on_actionD_viation_angles_triggered()
{
    deviationNormales(&mesh);
}



void MainWindow::on_actionBoite_englobante_triggered()
{
    QVector<float> result = boiteEnglobante(&mesh);
    QString text("Boîte englobante: ");
    for(int i = 0; i < result.size(); i++){
        text.append(QString::number(result[i]) + " ");
    }
    QMessageBox::information(this, tr("Boîte englobante"), text);
}

void MainWindow::on_actionBarycentre_de_la_forme_triggered()
{
    MyMesh::Point point = barycentreForme(&mesh);
    QMessageBox::information(this,
                             tr("Barycentre"),
                             QString("Le barycentre est de coordonnées (" + QString::number(point[0]) + ", " + QString::number(point[1]) + ", " + QString::number(point[2]) + ")"));
    this->displayMesh(&mesh, true);
}


void MainWindow::on_actionAngles_di_dres_triggered()
{
    std::vector<int> list_angle = anglesDihedres(&mesh);
    QString text("Angles dièdres : \n");
    for(int i = 0; i < list_angle.size(); i++){
        text.append(QString::number(i) + " : " + QString::number(list_angle[i]) + "\n");
    }
    QMessageBox::information(this, tr("Angles dièdres"), text);
}


void MainWindow::on_actionFr_quences_des_aires_triggered()
{
    std::vector<double> ratio = frequenceAires(&mesh);
    QString text("Fréquences des aires : \n");
    for(unsigned int i = 0; i < ratio.size()-1; i++){
        text.append("Aire entre " + QString::number(i*10) + "% et " + QString::number(i*10+10) + "% : " + QString::number(ratio[i]) + "% \n");
    }
    QMessageBox::information(this, tr("Fréquence aires"), text);

    ui->plainTextEdit->setPlainText("Aire totale du maillage = " + QString::number(ratio[ratio.size()-1]));
}

void MainWindow::on_actionFr_quences_valences_triggered()
{
    std::vector<double> ratio = frequenceVoisinageSommets(&mesh);
    QString text("Valences des sommets : \n");
    for(unsigned int i = 0; i < ratio.size()-1; i++){
        text.append("Valence de " + QString::number(ratio[ratio.size()-1] +i) + " : " + QString::number(ratio[i]) + "% \n");
    }
    QMessageBox::information(this, tr("Fréquence valences sommets"), text);
}

void MainWindow::on_actionV_rification_voisins_triggered()
{
    QVector<std::vector<int>> list = verificationVoisins(&mesh);

    QString text("Vérification voisinnage éléments : \n");

    if(list[0].size() == 0){
        text += "Toutes les faces ont au moins un voisin \n \n";
    }
    else{
        for(int i = 0; i<list[0].size(); i++){
            text += QString(list[0][i]);
        };
        text += "\n \n";
    }
    if(list[1].size() == 0){
        text += "Tous les points ont au moins une arête \n \n";
    }
    else{
        for(int i = 0; i<list[1].size(); i++){
            text += QString(list[1][i]);
        };
        text += "\n \n";
    }
    if(list[2].size() == 0){
        text += "Toutes les arêtes appartiennent à au moins une face";
    }
    else{
        for(int i = 0; i<list[2].size(); i++){
            text += QString(list[2][i]);
        };
        text += "\n \n";
    }

    QMessageBox::information(this, tr("Voisins"), text);
}

void MainWindow::on_actionCommandes_triggered()
{
    QMessageBox::information(this,
                             tr("Commandes"),
                             QString("- Tourner le maillage : clic gauche sur le maillage \n"
                                     "- Zoomer : molette de la souris OU clic gauche sur le maillage puis Ctrl OU clic gauche sur la maillage puis clic sur la molette \n"
                                     "- Translation dans le plan du viewer : clic gauche sur le maillage puis Alt"));
}
