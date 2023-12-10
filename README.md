# TextureProcessing

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
UV coordinates are distributed evenly along the boundary edges, achieved by calculating cumulative lengths of each segment of a boundary edge relative to its total length.

```

```
