#pragma once
#include <stdio.h>
#include "Cartesian3.h"
#include "AttributedObject.h"
#include "MatrixManager.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>


class DirectedEdgeMesh 
{

    public:
    //Coordinates of each vertice
    std::vector<Cartesian3> vertices;

    //color of vertices
    std::vector<Cartesian3> colours;

    //normals
    std::vector<Cartesian3> normals;
    
    // vector of texture coordinates
    std::vector<Cartesian3> textureCoords;

    // vector for the "first" directed edge for each vertex
    std::vector<long> firstDirectedEdge;

    //indices
    std::vector<unsigned int> faceVertices;

    // vector of the "other half" of the directed edges
    std::vector<long> otherHalf;

    //Boundaries of the mesh
    std::vector<unsigned int> boundaryEdge;

    //Ordered boundary edge
    std::vector<unsigned int> orderedBoundEdge;
    
    //Weights
    std::vector<std::vector<double>> weights;

    // centre of gravity - computed after reading
    Cartesian3 centreOfGravity;

    //Vertices that are actually used in mesh
    unsigned int validVertex;

    DirectedEdgeMesh();
	
    //Read vertices coordinates, normals, color, face indices from attributedObject
    //This AttributedObject has already read those data from an .obj file
	bool ReadFile(AttributedObject attriObject);

    //Texture parametrisation
    void Parameterize();

    //Get UV coordinate of start and end vertice
    std::vector<Cartesian3> GetUV(unsigned int currentEdgeIndex);


    //Calculate weights of every vertex
    void CalculateWeights();

    //Assign UV coordinates for boundary edge
    void BoundaryUvAssign();

    void InteriorUvAssign();

    void WriteFile();

    unsigned int getNext(unsigned int currentEdge);

};