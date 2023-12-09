#pragma once
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

class MatrixManager
{

public:

    const double epsilon = 1e-12;
    MatrixManager();

    vector<vector<double>> creatmatrix(int h,int l);

    vector<vector<double>> add(const vector<vector<double>>&A, const vector<vector<double>>&B);
    vector<vector<double>> minus(const vector<vector<double>>&A, const vector<vector<double>>&B);
    vector<vector<double>> multiply(const vector<vector<double>>&A, const vector<vector<double>>&B);

    vector<vector<double>> multiply_const(const vector<vector<double>>&A, double num);
    vector<vector<double>> trans(const vector<vector<double>>&A);
    vector<vector<double>> inverse(const vector<vector<double>>&A);

    void show(const vector<vector<double>>&A);
};