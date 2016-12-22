#include <iostream>
#include <fstream>
#include <stdio.h>
#include <Windows.h>
#include <conio.h>


using namespace std;
// e.g. f(x,y) = x^2 + y^2;
// define vector operations here(add Q*vect as well)


class vect
{
private:
		double *v;
		int dim;
public:
		inline vect()
		{
			v = NULL;
		}
		inline vect(int dimention)
		{
			
			v = new double[dimention];
			dim = dimention;
		}
		inline double &operator [] (unsigned i)
		{
			return v[i];
		}
		inline vect operator + (vect x)
		{
			vect res(dim);
			for(int i = 0; i < dim; i++)
			{
				res[i] = v[i] + x[i];
			}
			return res;
		}
		inline vect operator - (vect x)
		{
			vect res(dim);
			for(int i = 0; i < dim; i++)
			{
				res[i] = v[i] - x[i];
			}
			return res;
		}
		inline vect operator = (double *res)
		{
			for(int i = 0; i < dim; i++)
			{
				v[i] = res[i];
			}
			return *this;
		}

		inline vect operator * (double c)
		{
			vect res(dim);
			for(int i = 0; i < dim; i++)
			{
				res[i]=v[i] *c;
			}
			return res;
		}
		inline vect operator / (double c)
		{
			vect res(dim);
			for(int i = 0; i < dim; i++)
			{
				res[i]=v[i] /c;
			}
			return res;
		}
		inline double** operator * (vect b)//v*b; v - column; b - row;
		{
			double **res;
			res = new double*[dim];
			for(int i = 0; i < dim; i++)
			{
				res[i] = new double [dim];
				for(int j = 0; j < dim; j++)
				{
					res[i][j] = v[i]*b[j];
				}
			}
			return res;
		}

		double scal_mul(vect b)
		{
			double res=0;
			for (int i = 0; i < dim; i++)
			{
				res += v[i] * b[i];
			}
			return res;
		}
		double check_norm()
		{
			double res = 0;
			for (int i = 0; i < dim; i++)
			{
				res += pow(v[i],2);
			}
			return sqrt(res);
		}
};

class matrix
{
public:
	inline matrix(int dimention)
	{
		dim = dimention;
		M = new vect*[dim];
		for (int i = 0; i < dim; i++)
		{
			M[i] = new vect(dim);
		}
	}
	inline vect & matrix::operator [] (unsigned j)
	{
		return *M[j];
	}
	inline matrix operator = (double **X)
	{
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				
				*M[i] = X[i];
			}
		}
		return *this;
	}
	inline vect operator *(vect x)//A*b
	{
		double line_sum = 0;
		vect res(3);
		for (int i = 0; i < dim; i++)
		{
			line_sum = 0;
			for (int j = 0; j < dim; j++)
			{
				line_sum += (*M[i])[j] * x[j];	
			}
			res[i] = line_sum;
		}
		return res;
	}
	inline matrix operator *(double x)//A*b
	{
		matrix res(dim);
		for (int i = 0; i < dim; i++)
		{
			res[i] = *M[i] * x;
		}
		return res;
	}
	inline matrix operator +(matrix X)
	{
		matrix res(dim);
		for (int i = 0; i < dim; i++)
		{
				res[i] =*M[i] + X[i];
		}
		return res;
	}
	inline matrix operator -(matrix X)
	{
		matrix res(dim);
		for (int i = 0; i < dim; i++)
		{
			res[i] = *M[i] - X[i];
		}
		return res;
	}
	inline matrix operator /(double x)
	{
		matrix res(dim);
		for (int i = 0; i < dim; i++)
		{
			res[i] = *M[i] / x;
		}
		return res;
	}

	matrix transpose()
	{
		matrix res(dim);
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				res[i][j] = (*M[j])[i];
			}
		}
	return res;
	}
private:
	int dim;
	vect **M;
};
class QNM
{
public:
	QNM(char *filename)
	{
		//read data
		fstream fin;
		fin.open(filename);
		fin >> dim >> qnm_eps;
		//allocate
		Q = new double*[dim];
		x = new double [dim];

		for (int i = 0; i<dim; i++)
			fin >> x[i];
		fin.close();

		//Q is being initialised with Diag(dim) matrix
		for(int i=0; i<dim; i++)
		{
			Q[i] = new double [dim];
			for(int j=0; j<dim; j++) //i!=j ?? nvm
			{
				Q[i][j] = 0;
			}
			Q[i][i] = 1;
		}
	}
	~QNM()
	{
		delete x;
		for(int i = 0; i < dim; i++)
		{
			delete Q[i];
		}
		delete Q;
	}
	int check_diag_by_LLT(matrix Q_new);
	double Wulf(vect xk, vect dk);
	vect QNM_solution();
	double find_min(vect x,vect d, double a, double b);
private:
	int dim;
	double **Q; //Q0 matrix
	double *x;	//x0
	double qnm_eps;
	int max_iter;
};