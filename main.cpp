//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
//
//  September, 2020
//
//  -----------------------------
//  main.cpp
//  -----------------------------
//  
//  Loads assets, then passes them to the render window. This is very far
//  from the only way of doing it.
//  
////////////////////////////////////////////////////////////////////////

// system libraries
#include <iostream>
#include <fstream>

// QT
#include <QApplication>

// local includes
#include "RenderWindow.h"
#include "AttributedObject.h"
#include "DirectedEdgeMesh.h"
#include "RenderParameters.h"
#include "RenderController.h"

// main routine
int main(int argc, char **argv)
    { // main()
    // initialize QT
    QApplication renderApp(argc, argv);

    // check the args to make sure there's an input file
    if (argc != 2) 
        { // bad arg count
        // print an error message
        std::cout << "Usage: " << argv[0] << " geometry" << std::endl; 
        // and leave
        return 0;
        } // bad arg count

    //  use the argument to create a height field &c.
    AttributedObject AttributedObject;

    // open the input files for the geometry & texture
    std::ifstream geometryFile(argv[1]);

    // try reading it
    if (!(geometryFile.good()) || (!AttributedObject.ReadObjectStream(geometryFile)))
        { // object read failed 
        std::cout << "Read failed for object " << argv[1] << std::endl;
        return 0;
        } // object read failed

    // dump the file to out
//     AttributedObject.WriteObjectStream(std::cout);

    AttributedObject.Parameterize();
    

    // create some default render parameters
    RenderParameters renderParameters;

    // use the object & parameters to create a window
    RenderWindow renderWindow(&AttributedObject, &renderParameters, argv[1]);

    // create a controller for the window
    RenderController renderController(&AttributedObject, &renderParameters, &renderWindow);

    //  set the initial size
    renderWindow.resize(723, 580);

    // show the window
    renderWindow.show();

    // set QT running
    return renderApp.exec();

    } // main()
