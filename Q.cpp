#include "Q.h"

double f_x(vect x)
{
	//double result = 4 + pow(x[0] - 5.0, 4) + pow(x[1] - 6.0, 4);
	double result = pow(x[0] - 1, 2) + pow(x[1] - 1, 2);
	//double  result = sin(x[0] + x[1]) + pow(x[0] - x[1],2) - 1.5*x[0] + 2.5*x[1] + 1;
	//double result = x[0] + x[1]; // без глобального минимума
	//double result =  pow(x[0] - 1E-7, 2) + pow(x[1] - 1E-7, 2);
	return(result);
};

vect d_f(vect x, int dim) // equation is f(x,y) = ? Ask Лисицин.
{
	double result = 0; // init
	vect res(dim);
	for(int i = 0; i < dim; i++)
	{
		switch(i)
		{
		case 0:
			{
				//res[i] = 4*pow(x[i] - 5.0,3);
				res[i] = 2*(x[i] - 1);
				//res[i] = cos(x[0] + x[1]) + 2*(x[0] - x[1]) - 1.5;
				//res[i] = 1;
				//res[i] = 2*(x[i] - 1E-7);
				break;
			}
		case 1:
			{
				//res[i] = 4*pow(x[i] - 6.0,3);
				res[i] = 2*(x[i] - 1);
				//res[i] = cos(x[0] + x[1]) - 2*(x[0] - x[1]) + 2.5;
				//res[i] = 1;
				//res[i] = 2*(x[i] - 1E-7);
				break;
			}
		default:
			{
				res[i] = -1;
				printf("\nerror in d_f_x\n");
			}
		}
	}
	return res;
};

double QNM::Wulf(vect xk, vect dk)
{
	double eps1 = 10E-4, eps2 = 10E-3, alpha = 5, a_bot=0, a_top=0; //
	double thet1 = 2, thet2 = 0.5;
	vect xk_adk(dim);

	int stop_index = 0;
	while (stop_index < 1E4)
	{
		stop_index++;
		xk_adk = (dk*alpha) + xk;
		double temp1 = f_x(xk_adk) - f_x(xk) - eps1 * alpha * d_f(xk, dim).scal_mul(dk); // <0 true
		double temp2 = d_f(xk_adk, dim).scal_mul(dk) - eps2 * d_f(xk, dim).scal_mul(dk); // >0 true
		if (temp1 <= 0) // value1 - value2 < 0.0...01 ?
		{
			if (temp2 >= 0)
			{
				break;
			}
			else
			{
				a_bot = alpha;
				if (a_top == 0)
				{
					alpha *= thet1;
				}
			}
		}
		else
		{
			a_top = alpha;
			alpha = (1 - thet2)*a_bot + thet2*a_top;

		}
	}
	return alpha;
}

double QNM::find_min(vect x, vect dk,double a, double b) 
{
	const double eps = 1E-8;
        double e = 0;
        double x_1 = 0;
        double x_2 = 0;
		int i = 0;
        while (abs(b - a) >= eps) 
		{
			i++;
            e = (b - a) * 1E-5;
            x_1 = (b + a) / 2 - e;
            x_2 = (b + a) / 2 + e;

            if (f_x(x + dk * x_1) > f_x(x + dk *x_2  )) 
			{
                a = x_1;
            }
            else 
			{
                b = x_2;
            }
        }
		//printf("iter = %d ", i);
        return (a + b) / 2;
}

vect QNM::QNM_solution()
{
	fstream fout;
	fout.open("output.txt");

	vect d_k(dim), s_k(dim), r_k(dim), x_cur(dim), x_new(dim), rk_Qksk(dim);
	matrix Q_new(dim);
	matrix temp(dim);
	double alpha, scal_rk_sk = 0;
	max_iter = 0;
		x_cur = x;
		x_new = x_cur;

		Q_new = Q;
		 //x_k+1 = x_k + a*d_k
		int reset = 0;
		for(; max_iter < 1E+4 && d_f(x_new,dim).check_norm() >= qnm_eps ; x_cur = x_new, max_iter++, reset++)
		{
			d_k = d_f(x_cur, dim); // f'(x_k)
			d_k = Q_new * d_k;
			d_k = d_k*(-1); // d_k = -Q*f'
		//Wulf
			alpha = Wulf(x_cur,d_k);
			//alpha = this->find_min(x_cur,d_k,-10,10);
			x_new = d_k*alpha + x_cur;

			for(int i = 0; i < dim; i++)
				{
					printf("%f ",x_new[i]);
				}
				printf("\n");
		//Qk+1
			if (reset == dim || !check_diag_by_LLT(Q_new))
			{
				Q_new = Q;
				reset = -1;
			}
			else
			{
				r_k = x_new - x_cur;
				s_k = d_f(x_new, dim) - d_f(x_cur, dim);

				rk_Qksk = r_k - Q_new * s_k;
				scal_rk_sk = r_k.scal_mul(s_k);

				temp = rk_Qksk * r_k;
				Q_new = Q_new + (temp.transpose() + temp) / scal_rk_sk;
				temp = r_k *r_k;
				Q_new = Q_new - temp * rk_Qksk.scal_mul(s_k) / (scal_rk_sk*scal_rk_sk);
			}
		}

		printf("iterations = %d \n", max_iter);
		for(int i=0; i<dim; i++)
		{
			printf("%f ", x_cur[i]);
			fout << x_cur[i] <<" \n";
		}
	return x_cur;
}


int QNM::check_diag_by_LLT(matrix Q_new)
{
	double sum1, sum2;
	matrix LLT(dim);
	for (int i = 0; i < dim; i++)
	{
		sum1 = 0;
		for (int j = 0; j < i; j++)
		{
			sum2 = 0;
			for (int k = 0; k<j; k++)
			{
				sum2 += LLT[i][k] * LLT[j][k];
			}
			LLT[i][j] = (Q_new[i][j] - sum2) / LLT[j][j];
			sum1 += LLT[i][j] * LLT[i][j];
		}

		LLT[i][i] = sqrt(Q_new[i][i] - sum1);
		//check diag
		if (LLT[i][i] < 0)
			return 0;
	}

	return 1;
}