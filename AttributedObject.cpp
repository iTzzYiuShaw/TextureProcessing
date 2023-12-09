///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//	Variant on TexturedObject that stores explicit RGB
//	values for each vertex
//  
///////////////////////////////////////////////////

// include the header file
#include "AttributedObject.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

// include the Cartesian 3- vector class
#include "Cartesian3.h"

#define MAXIMUM_LINE_LENGTH 1024
#define REMAP_TO_UNIT_INTERVAL(x) (0.5 + (0.5*(x)))
#define REMAP_FROM_UNIT_INTERVAL(x) (-1.0 + (2.0*(x)))

#define N_ITERATIONS 100000
using namespace std;
// constructor will initialise to safe values
AttributedObject::AttributedObject()
    : centreOfGravity(0.0,0.0,0.0)
    { // AttributedObject()
    // force arrays to size 0
    vertices.resize(0);
    colours.resize(0);
    normals.resize(0);
    textureCoords.resize(0);
    firstDirectedEdge.resize(0);
    faceVertices.resize(0);
	otherHalf.resize(0);
    } // AttributedObject()

// read routine returns true on success, failure otherwise
bool AttributedObject::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
        // character to read
        char firstChar = geometryStream.get();
        
//         std::cout << "Read: " << firstChar << std::endl;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
            { // switch on first character
            case '#':       // comment line
                // read and discard the line
                geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
                break;
                
            case 'v':       // vertex data of some type
                { // some sort of vertex data
                // retrieve another character
                char secondChar = geometryStream.get();
                
                // bail if we ran out of file
                if (geometryStream.eof())
                    break;

                // now use the second character to choose branch
                switch (secondChar)
                    { // switch on second character
                    case ' ':       // space - indicates a vertex
                        { // vertex read
                        Cartesian3 vertex;
                        geometryStream >> vertex;
                        vertices.push_back(vertex);
//                         std::cout << "Vertex " << vertex << std::endl;
                        break;
                        } // vertex read
                    case 'c':       // c indicates colour
                        { // normal read
                        Cartesian3 colour;
                        geometryStream >> colour;
                        colours.push_back(colour);
//                         std::cout << "Colour " << colour << std::endl;
                        break;
                        } // normal read
                    case 'n':       // n indicates normal vector
                        { // normal read
                        Cartesian3 normal;
                        geometryStream >> normal;
                        normals.push_back(normal);
//                         std::cout << "Normal " << normal << std::endl;
                        break;
                        } // normal read
                    case 't':       // t indicates texture coords
                        { // tex coord
                        Cartesian3 texCoord;
                        geometryStream >> texCoord;
                        textureCoords.push_back(texCoord);
//                         std::cout << "Tex Coords " << texCoord << std::endl;
                        break;                  
                        } // tex coord
                    default:
                        break;
                    } // switch on second character 
                break;
                } // some sort of vertex data
                
            case 'f':       // face data
                { // face
				// make a hard assumption that we have a single triangle per line
                unsigned int vertexID;
                
                // read in three vertices
				for (unsigned int vertex = 0; vertex < 3; vertex++)
					{ // per vertex
					// read a vertex ID
					geometryStream >> vertexID;

					// subtract one and store them (OBJ uses 1-based numbering)
					faceVertices.push_back(vertexID-1);
					} // per vertex
				break;
                } // face
                
            // default processing: do nothing
            default:
                break;

            } // switch on first character

        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
                
            } // per vertex
        } // non-empty vertex set

// 	std::cout << "Centre of Gravity: " << centreOfGravity << std::endl;
// 	std::cout << "Object Size:       " << objectSize << std::endl;


    firstDirectedEdge.resize(vertices.size()+2,-1);
    otherHalf.resize(faceVertices.size()+2,-1);
    textureCoords.resize(vertices.size());
    normals.resize(vertices.size());
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


    // return a success code
    return true;
	} // ReadObjectStream()

// write routine
void AttributedObject::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
    geometryStream << "# " << (faceVertices.size()/3) << " triangles" << std::endl;
    geometryStream << std::endl;

    // output the vertex coordinates
    geometryStream << "# " << vertices.size() << " vertices" << std::endl;
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "v  " << std::fixed << vertices[vertex] << std::endl;

    // output the vertex colours
    geometryStream << "# " << colours.size() << " vertex colours" << std::endl;
    for (unsigned int vertex = 0; vertex < colours.size(); vertex++)
        geometryStream << "vc " << std::fixed << colours[vertex] << std::endl;

    // output the vertex normals
    geometryStream << "# " << normals.size() << " vertex normals" << std::endl;
    for (unsigned int vertex = 0; vertex < normals.size(); vertex++)
        geometryStream << "vn " << std::fixed << normals[vertex] << std::endl;

    // output the vertex coords
    geometryStream << "# " << textureCoords.size() << " vertex tex coords" << std::endl;
    for (unsigned int vertex = 0; vertex < textureCoords.size(); vertex++)
        geometryStream << "vt " << std::fixed << textureCoords[vertex] << std::endl;

    // and the faces
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        geometryStream << "f";
        
        // loop through # of vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
			{ // per vertex
            geometryStream << " ";
            geometryStream << faceVertices[face+vertex] + 1;
			} // per vertex
		// end the line
        geometryStream << std::endl;
        } // per face
    
    } // WriteObjectStream()

// routine to render
void AttributedObject::Render(RenderParameters *renderParameters)
    { // Render()
	// make sure that textures are disabled
	glDisable(GL_TEXTURE_2D);

	float scale = renderParameters->zoomScale;
	scale /= objectSize;
	// Scale defaults to the zoom setting
	glTranslatef(-centreOfGravity.x * scale, -centreOfGravity.y * scale, -centreOfGravity.z * scale);
		
	if (renderParameters->useWireframe)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	else
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    // start rendering
    glBegin(GL_TRIANGLES);
	
    // loop through the faces: note that they may not be triangles, which complicates life
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
        { // per face
        
		// now do a loop over three vertices
		for (unsigned int vertex = 0; vertex < 3; vertex++)
			{ // per vertex
			// set colour using vertex ID

            if(renderParameters->useTexCoords)
            {
                glColor3f
				(
				textureCoords[faceVertices[face+vertex]].x,
				textureCoords[faceVertices[face+vertex]].y,
				textureCoords[faceVertices[face+vertex]].z
				);

            }else if(renderParameters->renderNormalMap || renderParameters->useNormal)
            {
                //cout << normals[faceVertices[face+vertex]] << endl;
                glColor3f
				(
				normals[faceVertices[face+vertex]].x*1.5,
				normals[faceVertices[face+vertex]].y*1.5,
				normals[faceVertices[face+vertex]].z*1.5
				);
            }      
            else{

                glColor3f
				(
				colours[faceVertices[face+vertex]].x,
				colours[faceVertices[face+vertex]].y,
				colours[faceVertices[face+vertex]].z
				);
            }


			// use scaled xyz for vertex position
            if(renderParameters->renderTexture || renderParameters->renderNormalMap)
            {
                glVertex3f
				(
				scale * textureCoords[faceVertices[face+vertex]].x,
				scale * textureCoords[faceVertices[face+vertex]].y,
				scale * textureCoords[faceVertices[face+vertex]].z
				);
            }else
            {
                glVertex3f
				(
				scale * vertices[faceVertices[face+vertex]].x,
				scale * vertices[faceVertices[face+vertex]].y,
				scale * vertices[faceVertices[face+vertex]].z
				);
            }

			} // per vertex
        } // per face

    // close off the triangles
    glEnd();

    // revert render mode  
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    } // Render()


void AttributedObject::Parameterize()
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

    CalculateFaceNormal();

    CalculateVertexNormal();
}


void AttributedObject::InteriorUvAssign()
{

    unsigned delta = orderedBoundEdge.size(); 
    long iterate = 0;

    while (iterate < N_ITERATIONS)
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

}

unsigned int AttributedObject::getNext(unsigned int currentEdge)
{
    if(currentEdge % 3 == 2)
        return currentEdge-2;
    else
        return currentEdge+1;

}

std::vector<Cartesian3> AttributedObject::GetUV(unsigned int currentEdgeIndex)
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

void AttributedObject::CalculateWeights()
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

void AttributedObject::CalculateFaceNormal()
{
    unsigned int faceNum = faceVertices.size()/3;
    faceNormal.resize(faceNum);

    for(int i = 0; i < faceVertices.size();i++)
    {
        if(i % 3 == 2)
        {
            Cartesian3 v1 = vertices[faceVertices[i-2]];
            Cartesian3 v2 = vertices[faceVertices[i-1]];
            Cartesian3 v3 = vertices[faceVertices[i]];

            unsigned int faceID = i / 3;
            faceNormal[faceID] = ((v2-v1).cross(v3-v2)).unit();
        }
    }
}


void AttributedObject::CalculateVertexNormal()
{
    cout << "FaceNormal size: " << faceNormal.size() << endl;
    //Dealing with boundary vertex
    for(int i = 0; i < orderedBoundEdge.size(); i++)
    {
        unsigned int verID = faceVertices[orderedBoundEdge[i]];

        unsigned int firstEdge = firstDirectedEdge[verID];

        long curEdge = firstEdge;

        Cartesian3 normal;
        //Iterate all faces around this vertex
        //Left-hand side faces
        do
        {
            if(curEdge >= 0)
            {
                unsigned int faceID = curEdge / 3;
                normal = normal + faceNormal[faceID];

                curEdge = getNext(curEdge);
                curEdge = getNext(curEdge);
                curEdge = otherHalf[curEdge];
            }
            //cout << "CurEdge: " << curEdge << endl;

        }while(curEdge >= 0);

        //Right-hand side faces
        curEdge = firstEdge;
        do
        {
            unsigned int faceID = curEdge / 3;

            if(curEdge >= 0)
            {
                //Avoid adding the first faceNormal by twice
                if(curEdge != firstEdge)
                    normal = normal + faceNormal[faceID];
                
                curEdge = otherHalf[curEdge];

                if(curEdge >= 0)
                    curEdge = getNext(curEdge);
            }

        }while(curEdge >= 0);

        normal = normal.unit();
        normals[verID] = normal;
    }

    //Deal with internal vertices
    for(int i = 0; i < vertices.size(); i++)
    {
        //Make sure it's not a boundary vertex
        bool isBoundary = false;
        Cartesian3 normal;
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

        long firstEdge = firstDirectedEdge[i];

        if (firstEdge < 0)
            continue;
        
        unsigned int curEdge = firstEdge;
        do
        {
            unsigned int faceID = curEdge / 3;

            //cout <<"faceID: " <<faceID << endl;
            //cout <<"curEdge: " <<curEdge << endl;
            normal = normal + faceNormal[faceID];
            //cout << "face:" << normal << endl;

            curEdge = getNext(curEdge);
            curEdge = getNext(curEdge);
            curEdge = otherHalf[curEdge];

            
        } while (curEdge != firstEdge);
        
       
        normal = normal.unit();
        normals[i] = normal;
        //cout << "internal:" << normal << endl;
    }


}

void AttributedObject::BoundaryUvAssign()
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