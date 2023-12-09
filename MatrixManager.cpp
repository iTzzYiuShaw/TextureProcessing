#include "MatrixManager.h"
using namespace std;

MatrixManager::MatrixManager()
{}


vector<vector<double>> MatrixManager::creatmatrix(int h,int l)
{
	vector<vector<double>> v;
	for (int i = 0; i < h; i++)
	{
		vector<double>v1(l,0);
		v.push_back(v1);
	}
	return v;
}



vector<vector<double>> MatrixManager::add(const vector<vector<double>>&A, const vector<vector<double>>&B)
{
	int h=A.size();
	int l=A[0].size();

	vector<vector<double>> C;
	C=creatmatrix( h, l);
 
	for(int i=0;i<h;i++)
	{
		for (int j = 0; j < l; j++)
		{
			C[i][j]=A[i][j]+B[i][j];  
			if (abs(C[i][j])<epsilon)
			{
				C[i][j]=0.0;
			}
		}
	}
	return C;
}

vector<vector<double>> MatrixManager::minus(const vector<vector<double>>&A, const vector<vector<double>>&B)
{
	int h=A.size();
	int l=A[0].size();

	vector<vector<double>> C;
	C=creatmatrix( h, l);
 
	for(int i=0;i<h;i++)
	{
		for (int j = 0; j < l; j++)
		{
			C[i][j]=A[i][j]-B[i][j];  
			if (abs(C[i][j])<epsilon)
			{
				C[i][j]=0.0;
			}
		}
	}
	return C;
}

vector<vector<double>> MatrixManager::multiply(const vector<vector<double>>&A, const vector<vector<double>>&B)
{
	int A_h=A.size();
	int A_l=A[0].size();
	int B_h=B.size();
	int B_l=B[0].size();
	if(A_l !=B_h)
	{
		std::cout<<"These 2 matrices can not be muiltiplied with each other"<<std::endl;
		exit(0);
	}
	vector<vector<double>> C=creatmatrix(A_h,B_l);
	for (int i = 0; i < A_h; i++)
	{
		for (int j = 0; j < B_l; j++)
		{
			C[i][j]=0;
			for (int k = 0; k < A_l; k++)
			{
				C[i][j] +=A[i][k]*B[k][j];
			}
			if (abs(C[i][j])<epsilon)
			{
				C[i][j]=0.0;
			}
		}
	}
	return C;
}
vector<vector<double>> MatrixManager::multiply_const(const vector<vector<double>>&A, double num)
{
	int A_h=A.size();
	int A_l=A[0].size();
	vector<vector<double>> B=creatmatrix(A_h,A_l);
	for (int i = 0; i < A_h; i++)
	{
		for (int j = 0; j < A_l; j++)
		{
			B[i][j]=num*A[i][j];
		}
	}
	return B;
}

vector<vector<double>> MatrixManager::trans(const vector<vector<double>>&A)
{
	vector<vector<double>> AT=creatmatrix(A[0].size(),A.size());
	int h=AT.size();
	int l=AT[0].size();
	for (int i = 0; i <h ; i++)
	{
		for (int j = 0; j < l; j++)
		{
			AT[i][j]=A[j][i];
		}
	}
	return AT;
}

vector<vector<double>> MatrixManager::inverse(const vector<vector<double>>&A)
{
	if (A.size() != A[0].size())
	{
		std::cout<<"Illegal"<<std::endl;
		exit(0);
	}
	int n=A.size();
	vector<vector<double>> inv_A=creatmatrix(n,n);
	vector<vector<double>> L=creatmatrix(n,n);
	vector<vector<double>> U=creatmatrix(n,n);
	vector<vector<double>> inv_L=creatmatrix(n,n);
	vector<vector<double>> inv_U=creatmatrix(n,n);
	//LU factorisation

	for (int i = 0; i < n; i++)
	{
		L[i][i] = 1;   
	}
 
	for (int i = 0; i < n; i++)
	{
		U[0][i]=A[0][i];  
	}

	for (int i = 1; i < n; i++)
	{
		L[i][0]=1.0*A[i][0]/A[0][0];  
	}
 

	for (int i = 1; i < n; i++)
	{

		for (int j = i; j < n; j++)
		{
			double tem = 0;
			for (int k = 0; k < i; k++)
			{
				tem += L[i][k] * U[k][j];
			}
			U[i][j] = A[i][j] - tem;
			if (abs(U[i][j])<epsilon)
			{
				U[i][j]=0.0;
			}
		}
		
		for (int j = i ; j < n; j++)
		{
			double tem = 0;
			for (int k = 0; k < i; k++)
			{
				tem += L[j][k] * U[k][i];
			}
			L[j][i] = 1.0*(A[j][i] - tem) / U[i][i];
			if (abs(L[i][j])<epsilon)
			{
				L[i][j]=0.0;
			}
		}
 
	}

	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			if(i>j)
			{
				U[i][j]=0.0;
			}
			else if(i<j)
			{
				L[i][j]=0.0;
			}
		}
	}
 
	for (int i=0;i<n;i++) 
	{
		inv_U[i][i]=1/U[i][i];
		for (int k=i-1;k>=0;k--)
		{
			double s=0;
			for (int j=k+1;j<=i;j++)
			{
				s=s+U[k][j]*inv_U[j][i];
			}
			inv_U[k][i]=-s/U[k][k];
			if (abs(inv_U[k][i])<epsilon)
			{
				inv_U[k][i]=0.0;
			}
		}
	}

	for (int i=0;i<n;i++)  
	{
		inv_L[i][i]=1; 
		for (int k=i+1;k<n;k++)
		{
			for (int j=i;j<=k-1;j++)
			{
				inv_L[k][i]=inv_L[k][i]-L[k][j]*inv_L[j][i]; 
				if (abs(inv_L[k][i])<epsilon)
				{
					inv_L[k][i]=0.0;
				}
			}
		}
	}
	inv_A=multiply(inv_U,inv_L);
	return inv_A;
}

void MatrixManager::show(const vector<vector<double>>&A)
{
	int h=A.size();
	int l=A[0].size();
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < l; j++)
		{
			std::cout<<A[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
}