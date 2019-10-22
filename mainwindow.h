#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <QMessageBox>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include<string>

namespace Ui {
class MainWindow;
}

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

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    // TP1
    QVector<int> displayMeshStats(MyMesh* _mesh);
    QVector<std::vector<int>> verificationVoisins(MyMesh* _mesh);
    std::vector<double> frequenceAires(MyMesh* _mesh);
    std::vector<double> frequenceVoisinageSommets(MyMesh* _mesh);

    QVector<float> boiteEnglobante(MyMesh* _mesh);
    MyMesh::Point barycentreForme(MyMesh* _mesh);
    void showSelection(MyMesh* _mesh);
    void showFaceNormal(MyMesh* _mesh);
    void showVertexNormal(MyMesh* _mesh);
    void deviationNormales(MyMesh* _mesh);
    std::vector<int> anglesDihedres(MyMesh* _mesh);


    float faceArea(MyMesh* _mesh, int faceID);
    float baryArea(MyMesh* _mesh, int vertID);
    float angleFF(MyMesh *_mesh, int faceID0, int faceID1, int vertID0, int vertID1);
    float angleEE(MyMesh* _mesh, int vertexID, int faceID);
    void H_Curv(MyMesh* _mesh);
    void K_Curv(MyMesh* _mesh);

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

    void on_pushButton_h_clicked();

    void on_pushButton_k_clicked();

    void on_actionCharger_OBJ_triggered();

    void on_actionCourbure_Moyenne_triggered();

    void on_actionCourbure_Gaussienne_triggered();

    void on_actionD_viation_angles_triggered();

    void on_actionBoite_englobante_triggered();

    void on_actionBarycentre_de_la_forme_triggered();

    void on_actionAngles_di_dres_triggered();

    void on_actionFr_quences_des_aires_triggered();

    void on_actionFr_quences_valences_triggered();

    void on_actionV_rification_voisins_triggered();

    void on_actionCommandes_triggered();

private:

    bool modevoisinage;

    MyMesh mesh;
    std::string filename;

    QLabel *footerLabel;

    int vertexSelection;
    int edgeSelection;
    int faceSelection;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
