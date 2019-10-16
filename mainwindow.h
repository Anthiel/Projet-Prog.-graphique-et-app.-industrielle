#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <QMessageBox>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "util.h"

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // TP1
    void displayMeshStats(MyMesh* _mesh);
    void verificationVoisins(MyMesh* _mesh);
    void frequenceAires(MyMesh* _mesh);
    void frequenceVoisinageSommets(MyMesh* _mesh);

    QVector<float> boiteEnglobante(MyMesh* _mesh);
    MyMesh::Point barycentreForme(MyMesh* _mesh);



    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();

    void on_pushButton_bary_clicked();
    void on_pushButton_box_clicked();

private:

    bool modevoisinage;

    MyMesh mesh;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
