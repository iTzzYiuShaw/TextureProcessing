To compile:

module purge
module add legacy-eng
module load qt/5.15.2
qmake -project QT+=opengl
qmake
make

To run the program:

./TextureProcessing ../TextureProcessing/models/hamishbig.obj

In the models file, there is a .obj file called TestingFile1.obj that records all otherhalfs 
and boundary edges of target .obj file. Please check that to confirm.

On command line, the program will print out the number of boundary edges.
