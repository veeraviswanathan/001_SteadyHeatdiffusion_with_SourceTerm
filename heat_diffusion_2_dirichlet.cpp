#include<stdio.h>
#include<stdlib.h>
#include<array>
#include<cmath>


double L2norm(double *old, double *newarr);
int main()
{
	//declare domain variables
	int n_nodes = 100; // no of grid points
	double length = 1;// length of domain (m)
	double Tl = 600, Tr = 273; //left and right boundary condition (K)
	double Tinit = 200; // initial temperature for domain
	double k = 1; // Diffusivity (m2/s)
	double* T = nullptr; // Temperature array
	double* Told = nullptr; // Temperature array old
	double* S = nullptr; // Source array
	double residual = 1e20;
	double delta_x = length / n_nodes;
	double Fr, Fl;
	FILE *fp;

	fp=fopen("user_heat_diffusion_steady.out","w");
	fprintf(fp,"%s \t ","Niter");
	for (int kk=0;kk<n_nodes;kk++)
	{
		fprintf(fp,"%s[%d]\t","T",kk);
	}
	fprintf(fp,"%s \n","Residual");
	fclose(fp);


	//memory allocation
	T 		= (double*)calloc(n_nodes, sizeof(double));
	Told 	= (double*)calloc(n_nodes, sizeof(double));
	S 		= (double*)calloc(n_nodes, sizeof(double));

	// declare index variables
	int i;
	int iter = 0;
	int max_iter = 40;

	// initialisation block
	for (i = 0;i < n_nodes;i++)
	{
		T[i] = Tinit;
		Told[i] = T[i];
		S[i] = 0 * (delta_x * delta_x / k);
	}
	//S[int(n_nodes/2)]= 5000 * (delta_x * delta_x / k);
	//printf("%f %f \n", Tl, Tr);
	
	//Gauss Siedel iterative method
	fp=fopen("user_heat_diffusion_steady.out","a+");
	for (iter = 0;iter < max_iter;iter++)
	{

		// Assign old array
		for (i = 0;i < n_nodes;i++)
		{
			//printf("Temperature(%d) : %f \n",i,T[i]);
			Told[i] = T[i];
		}

		T[0] = (1.0/3)*(S[0] + (2 * Tl) + T[1]);
		for (i = 1;i < n_nodes - 2;i++)
		{
			T[i] = (1.0/2)*(S[i] + T[i + 1] + T[i - 1]);
		}
		T[n_nodes - 1] = (1.0/3)*(S[n_nodes - 1] + (2 * Tr) + T[n_nodes - 2]);

		//Flux at Right boundary
		Fr = (Tr - T[n_nodes-1]) / delta_x;
		//Flux at left boundary
		Fl = (Tl - T[0]) / delta_x;




		// L2norm error
		residual = L2norm(Told, T);
		//Write output file
		fprintf(fp,"%d \t",iter);
		for (i = 0;i < n_nodes;i++)
		{
			fprintf(fp,"%e \t",T[i]);
		}
		fprintf(fp,"%e \n",residual);
		printf("Iteration: %d Temperature midpoint: %e Residual: %f Flux Right: %e Flux left: %e\n", iter, T[int(n_nodes/2)], residual,Fr,Fl);
		if (residual < 1e-17)
		{
			printf("CONVERGED!");
			exit(0);
		}
		
	}
	fclose(fp);
	return(0);
}

double L2norm(double* Old, double* newarr)
{
	double len1 = sizeof(Old);
	double len2 = sizeof(newarr);
	double sum = 0.0;
	double dummy;
	if (len1 != len2)
	{
		printf("ERROR: previous and current arrays have different length - exiting");
			exit(0);
	}
	for (int i = 0;i < len1;i++)
	{
		dummy = (Old[i] - newarr[i]) * (Old[i] - newarr[i]);
		sum += sqrt(dummy);
	}
	return(sum);
}







