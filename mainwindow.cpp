#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>

/***************************** Constructor / Destructor *******************************/

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
    this->setWindowTitle("M2GIG - Projet PGAI");
}
MainWindow::~MainWindow()
{
    delete ui;
}

/***************************** TP1 *******************************/

// Affichage du nombre de sommets et de faces
void MainWindow::displayMeshStats(MyMesh* _mesh)
{
    qDebug() << "Ce maillage possède" << _mesh->n_vertices() << "sommets,"<< _mesh->n_faces() << "faces";
}
void MainWindow::verificationVoisins(MyMesh* _mesh)
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
void MainWindow::frequenceAires(MyMesh* _mesh)
{
    double maxA = Util::faceArea(&mesh, 0);
    double minA = Util::faceArea(&mesh, 0);
    double diff;
    std::vector <int> compt(10, 0);
    std::vector <double> ratio(10, 0);

    for(MyMesh::FaceIter curF = _mesh->faces_begin() ; curF != _mesh->faces_end() ; curF ++){
        if (Util::faceArea(&mesh, (*curF).idx()) > maxA){
            maxA = Util::faceArea(&mesh, (*curF).idx());
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) < minA){
            minA = Util::faceArea(&mesh, (*curF).idx());
        }
    }

    diff = maxA-minA;

    for(MyMesh::FaceIter curF = _mesh->faces_begin() ; curF != _mesh->faces_end() ; curF ++){
        if (Util::faceArea(&mesh, (*curF).idx()) > minA && Util::faceArea(&mesh, (*curF).idx()) <= minA+((10.0*diff)/100.0)){
            compt[0] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((10.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= minA+((20.0*diff)/100.0)){
            compt[1] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((20.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= minA+((30.0*diff)/100.0)){
            compt[2] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((30.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= minA+((40.0*diff)/100.0)){
            compt[3] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((40.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= minA+((50.0*diff)/100.0)){
            compt[4] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((50.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= minA+((60.0*diff)/100.0)){
            compt[5] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((60.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= minA+((70.0*diff)/100.0)){
            compt[6] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((70.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= minA+((80.0*diff)/100.0)){
            compt[7] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((80.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= minA+((90.0*diff)/100.0)){
            compt[8] ++;
        }
        else if (Util::faceArea(&mesh, (*curF).idx()) > minA+((90.0*diff)/100.0) && Util::faceArea(&mesh, (*curF).idx()) <= maxA){
            compt[9] ++;
        }
    }

    for (unsigned int i = 0; i < compt.size() ; i++){
        ratio[i] = (1.0*compt[i])/(1.0*_mesh->n_faces())*100;
        qDebug() << i*10 << "-" << i*10+10 << " aire : " << ratio[i] << "%";
    }
}
void MainWindow::frequenceVoisinageSommets(MyMesh* _mesh)
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

// Calcul moyenne des normales aux faces concourantes


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

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    this->resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    this->displayMesh(&mesh);

    this->displayMeshStats(&mesh);
    this->verificationVoisins(&mesh);
    this->frequenceAires(&mesh);
    this->frequenceVoisinageSommets(&mesh);
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

