
#include "DirectedEdgeMesh.h"
#include <iostream>
#include <fstream>
using namespace std;
DirectedEdgeMesh::DirectedEdgeMesh()
{

    vertices.resize(0);
    normals.resize(0);
    colours.resize(0);
    textureCoords.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    boundaryEdge.resize(0);
}


bool DirectedEdgeMesh::ReadFile(AttributedObject attriObject)
{

    //Deep copy of vectors
    vertices.resize(attriObject.vertices.size());
    normals.resize(attriObject.normals.size());
    colours.resize(attriObject.colours.size());
    faceVertices.resize(attriObject.faceVertices.size());
    textureCoords.resize(attriObject.vertices.size());

    vertices = attriObject.vertices;
    normals = attriObject.normals;
    colours = attriObject.colours;
    faceVertices = attriObject.faceVertices;

    firstDirectedEdge.resize(vertices.size()+2,-1);
    otherHalf.resize(faceVertices.size()+2,-1);

    //TODO: Looking for other half of each edge
    //TODO: Looking for firstDirectedEdge of each vertices
    //TODO: Looking for boundary edges

    int currentVertex = 0;
    int numOfUsedVer = 0;
    for(unsigned int i = 0; i < faceVertices.size(); i++)
    {
        
        if(i == faceVertices.size())
        {
            break;
        }

        //"To" vertice
        unsigned int to = faceVertices[i];
        unsigned int from;

        //Calculate the index of "From" vertex
        int delta = i % 3;
        if( delta == 0 )
            from = faceVertices[i + 2];
        else
            from = faceVertices[i - 1];

        otherHalf[i] = -2;

        /**
         * Goal 1:
         * For every edge that has already been checked
         * Go through every edge in every face that has already been checked
         * See if there is an other half edge
         * ---------------------------------------------------------
         * Goal 2:
         * Also, assign firstDirectedEdge for each vertex.
         * Notice: assign the first directed edge that STARTS FROM this vertex
        */
        for(int k = 0; k < i - delta; k++)
        {
            
            unsigned int F_to = faceVertices[k];
            unsigned int F_from;

            //Calculate the index of "From" vertex
            int F_delta = k % 3;
            if( F_delta == 0)
                F_from = faceVertices[k + 2];
            else
                F_from = faceVertices[k - 1];

            if(to == F_from && from == F_to)
            {
                //Assigns other half for both of the edges;
                //Notice: the index and value for otherHalf are both the indices of directed edges
                otherHalf[i] = k;
                otherHalf[k] = i;
            }
        }

        //Assigning first directed Edge
        if(firstDirectedEdge[from] == -1)
        {
            firstDirectedEdge[from] = i;
            numOfUsedVer++;
        }
    }


    cout << numOfUsedVer << endl;
    validVertex=numOfUsedVer;
    /**
     * If a directed edge has no otherhalf, it is a boundary edge
     * Loding the boundary edges.
    */
    for(unsigned int i = 0; i < otherHalf.size(); i++)
    {
        if(otherHalf[i] == -2)
        {
            boundaryEdge.push_back(i);
        }
    }

    /** Todo: Arrange the boundary, so that in boundaryEdge vector the boundaries are ordered
     * The order of boundary edge is consistent with the whole exterior boundary of the mesh. 
     * So that The boundary can be mapped on to the edges of the square
     * ------------------------------------------------
     * We will pick a randon boundary as a starting edge, iterate along with the exterior boundary of the mesh,
     * The "From" vertex of next boundary is the "To" vertex of this boundary
    */
    unsigned int cur_to_index = boundaryEdge[0];
    unsigned int start = cur_to_index;
    while(true)
    {    
        //Iterate every boundary, look for an edge whose "From" is "cur_to_index"
        for(int i = 0; i < boundaryEdge.size();i++)
        {

            unsigned int to_index = boundaryEdge[i];
            unsigned int from_index;

            //Get the index of "From" vertex
            int delta = to_index % 3;
            if( delta == 0 )
                from_index = to_index + 2;
            else
                from_index = to_index - 1;

            if(faceVertices[cur_to_index] == faceVertices[from_index])
            {
                //Push this edge into the array
                orderedBoundEdge.push_back(to_index);                
                //Travel to next boundary
                cur_to_index = to_index;               
                break;
            }
        }

        //Break if we travel back to the start edge
        if(cur_to_index == start)
            break;
    }
    std::cout << "Size of orderedBoundary: " << orderedBoundEdge.size() << std::endl;
    std::cout << "Size of Boundary: " << boundaryEdge.size() << std::endl;

    return true;
}


void DirectedEdgeMesh::Parameterize()
{

    BoundaryUvAssign();

    MatrixManager mm;
    
    //Fill the boundary data into matrix A. for every boundary edge, the weight of other point is 0
    //And the weight of itself is 1
    weights = mm.creatmatrix(vertices.size(),vertices.size());
    for(int i = 0; i < orderedBoundEdge.size(); i++)
    {
        //Index of Indices array
        unsigned int ind_Indices = orderedBoundEdge[i];

        //Index of Vertices array
        unsigned int ind_Vertices = faceVertices[ind_Indices];

        weights[i][i] = 1.0f;
    }

    //We start from where we ended in the last for loop
    //First directed edge: first directed edge that comes 'FROM' this vertex
    //Which gives us the first neibouring "To" vertex
    CalculateWeights();

    InteriorUvAssign();

}



void DirectedEdgeMesh::InteriorUvAssign()
{

    unsigned delta = orderedBoundEdge.size(); 
    long iterate = 0;
    long maxIteration = 10000000;

    while (iterate < maxIteration)
    {
        //Break the while loop if the difference is less than toleration
        double tolU = 0.0f;
        double tolV = 0.0f;
        unsigned int count = 0;
        for(int i = 0 ; i < vertices.size(); i++)
        {
            //Make sure it's not a boundary vertex
            bool isBoundary = false;
            for(int k = 0 ; k < orderedBoundEdge.size();k++)
            {
                unsigned int edgeIndex = orderedBoundEdge[k];
                unsigned int verIndex = faceVertices[edgeIndex];

                if (verIndex == i)
                {
                    isBoundary = true;
                    break;
                }          
            }
            if(isBoundary)
                continue;


            unsigned int startEdgeIndex = firstDirectedEdge[i];
            if(startEdgeIndex == -1)
                continue;

            Cartesian3 uv;
            Cartesian3 formerUV = textureCoords[i];
            unsigned int currentEdge_To = startEdgeIndex;
            do
            {
                unsigned int jIndex = faceVertices[currentEdge_To];
                uv = uv + (weights[delta + count][jIndex] * textureCoords[jIndex]) ;

                //Move to next neighbouring vertex
                currentEdge_To = getNext(currentEdge_To);
                currentEdge_To = getNext(currentEdge_To);
                currentEdge_To = otherHalf[currentEdge_To];

            } while (currentEdge_To != startEdgeIndex);
            
            tolU += abs(uv.x - formerUV.x);
            tolV += abs(uv.y - formerUV.y);

            textureCoords[i] = uv;
            count++;
        }

        iterate++;
        tolU /= count;
        tolV /= count;

        if (tolU < 1e-5 && tolV < 1e-5)
        {
            break;
        }  
    }

    // for(int i = 0 ; i < textureCoords.size(); i++)
    // {

    //     if(firstDirectedEdge[i] == -1)
    //         continue;

    //     bool isBoundary = false;
    //     for(int k = 0 ; k < orderedBoundEdge.size();k++)
    //     {
    //         unsigned int edgeIndex = orderedBoundEdge[k];
    //         unsigned int verIndex = indices[edgeIndex];

    //         if (verIndex == i)
    //         {
    //             isBoundary = true;
    //             break;
    //         }          
    //     }
        
    //     if(!isBoundary)
    //         cout << i << " : " << textureCoords[i].x << "  " << textureCoords[i].y << endl;
    // }

}

unsigned int DirectedEdgeMesh::getNext(unsigned int currentEdge)
{
    if(currentEdge % 3 == 2)
        return currentEdge-2;
    else
        return currentEdge+1;

}
std::vector<Cartesian3> DirectedEdgeMesh::GetUV(unsigned int currentEdgeIndex)
{
    Cartesian3 UV_start;
    Cartesian3 UV_end;

    if(currentEdgeIndex==0)
    {
        UV_start.x = 0.0f;
        UV_start.y = 0.0f;

        UV_end.x = 0.0f;
        UV_end.y = 1.0f;
    }else if(currentEdgeIndex == 1)
    {
        UV_start.x = 0.0f;
        UV_start.y = 1.0f;

        UV_end.x = 1.0f;
        UV_end.y = 1.0f;
    }else if(currentEdgeIndex == 2)
    {
        UV_start.x = 1.0f;
        UV_start.y = 1.0f;

        UV_end.x = 1.0f;
        UV_end.y = 0.0f;
    }else if(currentEdgeIndex == 3)
    {
        UV_start.x = 1.0f;
        UV_start.y = 0.0f;

        UV_end.x = 0.0f;
        UV_end.y = 0.0f;
    }
    std::vector<Cartesian3> UV;
    UV.push_back(UV_start);
    UV.push_back(UV_end);

    return UV;         
}

void DirectedEdgeMesh::CalculateWeights()
{
    unsigned int start = orderedBoundEdge.size();
    unsigned int count = 0;
    std::vector<unsigned int> invalidVerIndx;
    for(int i = 0; i < vertices.size(); i++)
    {
        //Make sure it's not a boundary vertex
        bool isBoundary = false;
        for(int k = 0 ; k < orderedBoundEdge.size();k++)
        {
            unsigned int edgeIndex = orderedBoundEdge[k];
            unsigned int verIndex = faceVertices[edgeIndex];

            if (verIndex == i)
            {
                isBoundary = true;
                break;
            }          
        }
        if(isBoundary)
            continue;

        //To get the next edge (other half) (Index of the Edge)
        unsigned int startEdgeIndex = firstDirectedEdge[i];
        if(startEdgeIndex == -1)
        {
            invalidVerIndx.push_back(i);
            continue;
        }
            
        
        //Iterate all neibouring vertices
        unsigned int currentEdge_To = startEdgeIndex;
        double accumWeight = 0.0f;

        //Store index of neighbouring vertices
        std::vector<unsigned int> targetIndex;
        do
        {
            //Get coordinate of vertex that is in the middle
            unsigned int jIndex = faceVertices[currentEdge_To];

            Cartesian3 vj = vertices[jIndex];
            Cartesian3 vi = vertices[i];

            //Todo: Get next neighbouring vertex
            unsigned nextEdge_to;

            if(currentEdge_To % 3 == 2)
                nextEdge_to = currentEdge_To-2;
            else
                nextEdge_to = currentEdge_To+1;

            unsigned nextVer = faceVertices[nextEdge_to];
            Cartesian3 nextVerCoord = vertices[nextVer];

            //Todo: Get former neighbouring vertex
            //1:OtherHalf of currentVer
            //2:next

            unsigned otherHalf_cur_to = otherHalf[currentEdge_To];
            unsigned F_nextEdge_to;
            if(otherHalf_cur_to % 3 == 2)
                F_nextEdge_to = otherHalf_cur_to-2;
            else
                F_nextEdge_to = otherHalf_cur_to+1;

            unsigned formerVer = faceVertices[F_nextEdge_to];
            Cartesian3 formerVerCoord = vertices[formerVer];

            //Todo: Calculate weight
            double dis_ij = (vi-vj).length();

            //tan(&_ij/2)
            Cartesian3 b1 = vj-vi;
            Cartesian3 a1 = formerVerCoord-vi;
            double cos1 = (a1.dot(b1)) / (a1.length() * b1.length());
            double sin1 = sqrt(1.0 - pow(cos1,2.0));
            double tan1 = sin1 / (1.0 + cos1);

            //tan(y_ij/2)
            Cartesian3 b2 = vj-vi;
            Cartesian3 a2 = nextVerCoord - vi;
            double cos2 = (a2.dot(b2)) / (a2.length() * b2.length());
            double sin2 = sqrt(1.0 - pow(cos2,2.0));
            double tan2 = sin2 / (1.0 + cos2);

            double weight_ij = (tan1 + tan2) / dis_ij;
            
            //Todo: Store this weight value into our matrix A
            //1: Index of vertex i in matrix A
            //2: Index of vertex j

            unsigned int i_index_A = start+count;
            weights[i_index_A][jIndex] = weight_ij;
            accumWeight+=weight_ij;

            targetIndex.push_back(jIndex);

            //Move to next neighbouring vertex: otherhalf(next->next)
            if(nextEdge_to % 3 == 2)
                nextEdge_to = nextEdge_to-2;
            else
                nextEdge_to = nextEdge_to+1;
            
            currentEdge_To = otherHalf[nextEdge_to];

            if(currentEdge_To == startEdgeIndex)
            {

                for(int j = 0 ; j < targetIndex.size();j++)
                {
                    unsigned index = targetIndex[j];
                    weights[i_index_A][index] = weights[i_index_A][index]/accumWeight;
                    //cout << " " <<weights[i_index_A][index];
                }
                //cout << " " << endl;
            }

        } while (currentEdge_To != startEdgeIndex);

        count++;
    }
    
}

void DirectedEdgeMesh::WriteFile()
{
    ofstream objFile;
    objFile.open("../TextureProcessing/models/hamish3k.obj", ios::out);

    if (objFile.is_open())
    {
        objFile << "# " << (faceVertices.size()/3) << " triangles" << std::endl;
        objFile << std::endl;

                objFile << "# " << vertices.size() << " vertices" << std::endl;
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            objFile << "v  " << std::fixed << vertices[vertex] << std::endl;

        // // output the vertex colours
        // objFile << "# " << colours.size() << " vertex colours" << std::endl;
        // for (unsigned int vertex = 0; vertex < colours.size(); vertex++)
        //     objFile << "vc " << std::fixed << colours[vertex] << std::endl;

        // // output the vertex normals
        // objFile << "# " << normals.size() << " vertex normals" << std::endl;
        // for (unsigned int vertex = 0; vertex < normals.size(); vertex++)
        //     objFile << "vn " << std::fixed << normals[vertex] << std::endl;

        // // output the vertex coords
        // objFile << "# " << textureCoords.size() << " vertex tex coords" << std::endl;
        // for (unsigned int vertex = 0; vertex < textureCoords.size(); vertex++)
        //     objFile << "vt " << std::fixed << textureCoords[vertex] << std::endl;

        // and the faces
        for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
            objFile << "f";
            
            // loop through # of vertices
            for (unsigned int vertex = 0; vertex < 3; vertex++)
                { // per vertex
                objFile << " ";
                objFile << faceVertices[face+vertex] + 1;
                } // per vertex
            // end the line
            objFile << std::endl;
        } // per face
        
        //Other Half
        for(unsigned int i = 0; i < otherHalf.size();i++)
        {
            objFile << "otherHalf";
            objFile << " ";
            objFile << i;

            objFile << "  ";
            objFile << otherHalf[i];

            objFile << std::endl;
        }
        objFile << std::endl;

        for(unsigned int i = 0; i < firstDirectedEdge.size();i++)
        {
            objFile << "FirstDirectedEdge";
            objFile << " ";
            objFile << i;

            objFile << "  ";
            objFile << firstDirectedEdge[i];

            objFile << std::endl;
        }

    }
    objFile.close();
} // WriteObjectStream()

void DirectedEdgeMesh::BoundaryUvAssign()
{
    //Todo: Assign an uv coordinate for every boundary edges
    unsigned int numOnEdge = (orderedBoundEdge.size() + 4) / 4;

    unsigned int n = (orderedBoundEdge.size() + 4) % 4;
    unsigned int m = 4 - n;

    unsigned int currentEdgeIndex = 0;

    // //Dealing with edges that have ((n+4)/4)+1 vertices
   for(int i = 0; i < n; i++)
    {
        unsigned int startIndex = currentEdgeIndex * numOnEdge;
        unsigned int endIndex = startIndex + numOnEdge;

        //Calculate total length from start edge to end edge
        double length = 0.0f;
        
        for(int i = startIndex; i < endIndex; i++)
        {
            //Thoes are indices of indices
            unsigned int to_index = orderedBoundEdge[i];
            unsigned int from_index;

            int delta = to_index % 3;
            if( delta == 0 )
                from_index = to_index + 2;
            else
                from_index = to_index - 1;

            unsigned int vToIndex = faceVertices[to_index];
            unsigned int vFromIndex = faceVertices[from_index];

            double cur_length = (vertices[vToIndex]-vertices[vFromIndex]).length();
            length+=cur_length;
        }

        //Get UV coordinates of 2 ends of a square's edge
        std::vector<Cartesian3> UV = GetUV(currentEdgeIndex);
        Cartesian3 UV_start = UV[0];
        Cartesian3 UV_end = UV[1];

        double tempLength = 0;
        for(int i = startIndex; i < endIndex; i++)
        {

            unsigned int to_index = orderedBoundEdge[i];
            unsigned int from_index;

            int delta = to_index % 3;
            if( delta == 0 )
                from_index = to_index + 2;
            else
                from_index = to_index - 1;

            unsigned int vToIndex = faceVertices[to_index];
            unsigned int vFromIndex = faceVertices[from_index];

            if (i == startIndex || i == endIndex-1)
            {
                
                textureCoords[vToIndex] = UV_start;
                if(i == endIndex-1)
                {
                    textureCoords[vToIndex] = UV_end;             
                }

                //std::cout << "UV: " << textureCoords[vToIndex] << std::endl;
                continue;
            }
            
            double cur_length = (vertices[vToIndex]-vertices[vFromIndex]).length();
            tempLength+=cur_length;

            Cartesian3 curUV = (tempLength / length) * (UV_end-UV_start) + UV_start;
            textureCoords[vToIndex] = curUV;

            //std::cout << "UV: " << textureCoords[vToIndex] << std::endl;
        }

        // std::cout << "startIndex: "<< startIndex << "endIndex: " << endIndex << std::endl;
        // std::cout << "Total length: " << length << std::endl;
        // std::cout << "Num of this edge: " << endIndex - startIndex + 1 << std::endl;
        // std::cout << "----------------------------------------------" << std::endl;

        currentEdgeIndex++;
    }

    //Dealing with edges that have ((n+4)/4) vertices
    for(int i = 0; i < m; i++)
    {
        unsigned int startIndex = n * (numOnEdge) + i * (numOnEdge-1);
        unsigned int endIndex = startIndex + numOnEdge-1;

        //Calculate total length from start edge to end edge
        double length = 0.0f;
        
        for(int i = startIndex; i < endIndex; i++)
        {
            //Thoes are indices of indices
            unsigned int to_index = orderedBoundEdge[i];
            unsigned int from_index;

            int delta = to_index % 3;
            if( delta == 0 )
                from_index = to_index + 2;
            else
                from_index = to_index - 1;

            unsigned int vToIndex = faceVertices[to_index];
            unsigned int vFromIndex = faceVertices[from_index];

            double cur_length = (vertices[vToIndex]-vertices[vFromIndex]).length();
            length+=cur_length;
        }


        //Get UV coordinates of 2 ends of a square's edge
        std::vector<Cartesian3> UV = GetUV(currentEdgeIndex);
        Cartesian3 UV_start = UV[0];
        Cartesian3 UV_end = UV[1];

        double tempLength = 0;
        for(int i = startIndex; i < endIndex; i++)
        {

            unsigned int to_index = orderedBoundEdge[i];
            unsigned int from_index;

            int delta = to_index % 3;
            if( delta == 0 )
                from_index = to_index + 2;
            else
                from_index = to_index - 1;

            unsigned int vToIndex = faceVertices[to_index];
            unsigned int vFromIndex = faceVertices[from_index];

            if (i == startIndex || i == endIndex-1)
            {
                
                textureCoords[vToIndex] = UV_start;
                if(i == endIndex-1)
                {
                    textureCoords[vToIndex] = UV_end;             
                }

                //std::cout << "UV: " << textureCoords[vToIndex] << std::endl;
                continue;
            }

            double cur_length = (vertices[vToIndex]-vertices[vFromIndex]).length();
            tempLength+=cur_length;

            Cartesian3 curUV = (tempLength / length) * (UV_end-UV_start) + UV_start;
            textureCoords[vToIndex] = curUV;

            //std::cout << "UV: " << textureCoords[vToIndex] << std::endl;
        }

        //std::cout << "startIndex: "<< startIndex << "endIndex: " << endIndex << std::endl;
        //std::cout << "Total length: " << length << std::endl;
        //std::cout << "Num of this edge: " << endIndex - startIndex + 1 << std::endl;
        //std::cout << "----------------------------------------------" << std::endl;

        currentEdgeIndex++;
    }
    //End: Assign an uv coordinate for every boundary edges
}