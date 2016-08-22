#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CUTE _stdcall

int CUTE pca(int n, int T, double *pricedata, double *x, double *Q)
{
	double *mean, *covariance, *returns, *w0, *w1, *w, *V, *F, *D;
	int retcode = 1;
	mean = (double *)calloc(n, sizeof(double));
	returns = (double *)calloc(n*(T - 1), sizeof(double));
	covariance = (double *)calloc(n*n, sizeof(double));
	w = (double *)calloc(n, sizeof(double));
	w0 = (double *)calloc(n, sizeof(double));
	w1 = (double *)calloc(n, sizeof(double));
	V = (double *)calloc(n * 10, sizeof(double));
	F = (double *)calloc(10 * 10, sizeof(double));
	D = (double *)calloc(n*n, sizeof(double));
	//retmat = (double *)calloc(10 * 10 + n * 10 + n*n, sizeof(double));
	double sum, sum1, norm, b, c, lambda, wTw0;

	int i, j, k;

	//Get Returns From Price Data
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < T - 1; j++)
		{
			returns[i*(T - 1) + j] = (pricedata[i*T + j + 1] - pricedata[i*T + j]) / pricedata[i*T + j];
		}
	}

	//Get Mean Return for each asset
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < T - 1; j++)
		{
			mean[i] += returns[i*(T - 1) + j] / (T - 1);
		}
	}

	//Get Covariance matrix for returns data
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < T - 1; k++)
			{
				covariance[i*n + j] += ((returns[i*(T - 1) + k] - mean[i]) * (returns[j*(T - 1) + k] - mean[j]))/(T-2);
				Q[i*n + j] = covariance[i*n + j];
			}
		}
	}

	//Generate Random Vector
	sum = 0; 
	for (i = 0; i < n; i++)
	{
		w[i] = rand() % 2;
		sum += pow(w[i], 2); 
	}
	norm = pow(sum, .5);

	
	for (k = 0; k < 10; k++)
	{
		c = 10000;
		for (i = 0; i < n; i++)
		{
			w1[i] = w[i];
		}
		//Power Method for Obtaining Eigenvector
		while (c > .00000001)
		{
			for (i = 0; i < n; i++)
			{
				w0[i] = w1[i];
				w1[i] = 0;
			}
			
			sum = 0;
			for (i = 0; i < n; i++)
			{ 
				for (j = 0; j < n; j++)
				{
					w1[i] += covariance[i*n + j] * w0[j];
				}
				sum += pow(w1[i], 2);
			}
			norm = pow(sum, .5);

			for (i = 0; i < n; i++)
			{
				w1[i] = w1[i] / norm;
			}
			b = c;
			c = 0;
			for (i = 0; i < n; i++)
			{
				c += fabs(w1[i] - w0[i]);
			}

		}
		
		sum1 = 0;
		for (i = 0; i < n; i++)
		{
			sum1 += pow(w1[i], 2);
		}
		norm = pow(sum1, .5);
		sum = 0;
		for (i = 0; i < n; i++)
		{
			V[i * 10 + k] = w1[i]/norm;
			sum += (w1[i]/norm) * covariance[i];
		}
		lambda = sum / w1[0];
		F[k * 10 + k] = lambda;
		wTw0 = 0;
		for (i = 0; i < n; i++)
		{
			wTw0 += w1[i] * w[i];
		}
		for (i = 0; i < n; i++)
		{
			w[i] = w[i] - wTw0*w1[i];
		}

		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				covariance[i*n + j] = covariance[i*n + j] - lambda*w1[i] * w1[j];
			}
		}

	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < 10; k++)
			{
				D[i*n + j] += V[i * 10 + k] * V[j * 10 + k] * F[k * 10 + k];
			}
			D[i*n + j] = covariance[i*n + j] - D[i*n + j];
		}
	}

	for (i = 0; i < 10; i++)
	{
		for (j = 0; j < 10; j++)
		{
			x[i * 10 + j] = F[i * 10 + j];
		}
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 10; j++)
		{
			x[10 * 10 + i * 10 + j] = V[i * 10 + j];
		}
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			x[10 * 10 + 10 * n + i*n + j] = D[i*n + j];
		}
	}



	free(returns);
	free(mean);
	free(covariance);
	free(w0);
	free(w1);
	free(V);
	free(F);
	free(D);

	return retcode;
}