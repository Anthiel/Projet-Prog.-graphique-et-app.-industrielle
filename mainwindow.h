#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <QMessageBox>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include<string>

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
    std::vector<double> frequenceAires(MyMesh* _mesh);
    std::vector<double> frequenceVoisinageSommets(MyMesh* _mesh);

    QVector<float> boiteEnglobante(MyMesh* _mesh);
    MyMesh::Point barycentreForme(MyMesh* _mesh);
    void showSelection(MyMesh* _mesh);
    void showFaceNormal(MyMesh* _mesh);
    void showVertexNormal(MyMesh* _mesh);
    void deviationNormales(MyMesh* _mesh);
    void anglesDihedres(MyMesh* _mesh);

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();

    void on_pushButton_bary_clicked();
    void on_pushButton_box_clicked();

    void on_vertexSelect_valueChanged(int arg1);

    void on_edgeSelect_valueChanged(int arg1);

    void on_faceSelect_valueChanged(int arg1);
    
    void on_pushButton_dev_clicked();

    void on_pushButton_dihedral_clicked();

    void on_pushButton_aire_clicked();

    void on_pushButton_freq_valence_clicked();

private:

    bool modevoisinage;

    MyMesh mesh;
    std::string filename;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
