# TextureProcessing


![efe780e829010dc25d9baab7eacf6ef](https://github.com/iTzzYiuShaw/TextureProcessing/assets/110170509/3fd23f5f-c52f-4066-99c7-34782b4970b2)

![75f961be6ea129ced04913240e0c95c](https://github.com/iTzzYiuShaw/TextureProcessing/assets/110170509/c9d1d9c4-91c7-4ec0-b5e8-8ff1a1f16729)


## Data pre-processing

### Looking for other-half of each edge
```
  // Iterate over all vertices that form faces

  // Calculate "From vertex" of that directed edge (The value of the directed edge is the id of STARTING VERTEX)

  //Go through every edge in every face that has already been checked, See if there is an other-half edge

  //Also, assign firstDirectedEdge for each vertex.
```

### Looking for boundary edges

```
//If a directed edge has no other half, it is a boundary edge
```

### Arrange the boundary, so that in boundary edge vector the boundaries are ordered
```
//We will pick a random boundary as a starting edge, iterate along with the exterior boundary of the mesh,
//The "From" vertex of the next boundary is the "To" vertex of this boundary
```

## Assigning boundary uv

For each boundary edge, the function maps the actual lengths (distances in physical space) between vertices to UV coordinates. 
This is done through linear interpolation, with start and end points provided by the GetUV function.
UV coordinates are distributed evenly along the boundary edges, achieved by calculating the cumulative lengths of each segment of a boundary edge relative to its total length.

```
    //Todo: Assign an uv coordinate for every boundary edges
    unsigned int numOnEdge = (orderedBoundEdge.size() + 4) / 4;

    unsigned int n = (orderedBoundEdge.size() + 4) % 4;
    unsigned int m = 4 - n;

    unsigned int currentEdgeIndex = 0;

// The number of vertices on each edge of the texture should be as even as possible,
// Therefore, the number of vertices on each edge could be either (int)(orderedBoundEdge.size() + 4) / 4 or (int)(orderedBoundEdge.size() + 4) / 4 + 1
// For example: if 14 vertices form the boundary of a mesh, then there should be (14+4)/4 = 4.5 vertices on each edge, which should be rounded into 4. Then the number of edges that own
// 4 vertices is 18 % 4 = 2. Then 2 edges own 4 vertices, and 2 edges own 5 vertices. To verify, there are going to be 2 * 4 + 2 * 5 = 18 vertices in total

            double cur_length = (vertices[vToIndex]-vertices[vFromIndex]).length();
            tempLength+=cur_length;

            Cartesian3 curUV = (tempLength / length) * (UV_end-UV_start) + UV_start;
            textureCoords[vToIndex] = curUV;
//achieved by calculating the cumulative lengths of each segment of a boundary edge relative to its total length.
```
## Calculate Weights Matrix

The function iterates over all neighboring vertices of a given vertex. This is done by cycling through connected edges.
For each neighboring vertex, it calculates the weight between the current vertex (vi) and its neighbor (vj).

## InteriorUvAssign

Iterative Approach:

The function uses an iterative approach to gradually converge towards the solution. The number of iterations is controlled by N_ITERATIONS, and the process is repeated until a convergence criterion is met or the maximum number of iterations is reached.
Convergence Criterion:

The convergence is checked using the variables tolU and tolV, which measure the change in UV coordinates between iterations. If the change falls below a certain threshold (1e-5 in this case), it is assumed that the solution has converged, and the loop breaks.
Updating Non-Boundary Vertices:

The main loop iterates over all vertices, skipping boundary vertices (determined by checking orderedBoundEdge and face vertices). For each non-boundary vertex, the function calculates new UV coordinates based on the weighted average of its neighboring vertices' UV coordinates.
The weights used here are presumably the same as those calculated in your previous function (CalculateWeights).
Calculation of New UV Coordinates:

For each non-boundary vertex, the function calculates a new UV coordinate (uv) as the weighted sum of the UV coordinates of its neighboring vertices. The weights are obtained from the weights matrix.
This process effectively smoothes out the UV mapping by averaging the influence of each vertex's neighbors, leading to a more evenly distributed UV map.
Tolerance Calculation:

The change in UV coordinates (tolU and tolV) is calculated as the sum of absolute differences between the new and former UV coordinates for all non-boundary vertices.
These values are averaged over the number of non-boundary vertices and used to check if the solution has converged.
