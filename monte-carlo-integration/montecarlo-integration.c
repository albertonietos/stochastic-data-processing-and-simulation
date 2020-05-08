/* Ordinary Monte Carlo integration
By Alberto Nieto
2018-10-19
*/
#include <math.h> //for mathfunctions
#include <stdio.h>
#include <stdlib.h>
#include <time.h> // for timing


double ordinary_montecarlo(double n)
    {
    	//seed rand function
        srand(time(NULL));
        // Define parameters to be used in the code
        double x=0, y=0, S_n=0, errorterm=0, error=0;
	    // loop according to sample size n
    	for (int i = 0; i < n; i++) {
	    		x = rand() / (RAND_MAX + 1.0); //samples
		    	y = 4/(1+x*x); // approximation

			    S_n= S_n+y; // uppdate for each loop according to sample size
			    errorterm= y*y + errorterm;
    	}
	    S_n= S_n/n; // Divide by sample size to get the approximated value
	    error=sqrt(errorterm/n-(S_n*S_n))/sqrt(n); // approximated error
    printf("Ordinary Monte Carlo Integration:\n");
    printf("S_n= %f, Estimated error = +/- %f, Real error = %f\n", S_n, error, (M_PI-S_n));
    }


double importance_sampling(double n)
	{
		//seed rand function
	 srand48(time(NULL)); //Make a new seed for the random number generator
    double x=0, y=0, S_n=0, errorterm=0, error=0, P=0, z=0;
		// loop according to sample size n
		for (int i = 0; i < n; i++) {

				z = drand48(); // randomize
				x = 2-sqrt(4-3*z); // define sample
				P = (4 - 2*x)/3; // define the p(x) function
				y = (4/(1+x*x)); // define the y function

				S_n += y/P; // do opertion and update according to sample size
				errorterm += ((y/P)*(y/P));
		}
		S_n = S_n/n;
		error = sqrt(errorterm/n - S_n*S_n)/sqrt(n);
	printf("Importance Sampling:\n");
	printf("S_n= %f, Estimated error +/- %f, Real error =%f\n", S_n, error, (M_PI-S_n));
	}


double control_variates(double n)
		{
			//seed rand function
		 srand48(time(NULL)); //Make a new seed for the random number generator
		     double x=0, y=0, S_n=0, errorterm=0, error=0, g=0, I=0;
			// loop according to sample size n
			for (int i = 0; i < n; i++) {

					x = drand48();
					//x = 2-sqrt(4-3*z);
					g = (4 - 2*x);
					y = (4/(1+x*x));
					I = 3.0;

					S_n += y-g+I;
					errorterm += ((y-g+I)*(y-g+I));
			}
			S_n = S_n/n;
			error = sqrt(errorterm/n - S_n*S_n)/sqrt(n);
		printf("Control Variates:\n");
		printf("S_n= %f, Estimated error +/- %f, Real error =%f\n", S_n, error, (M_PI-S_n));
		}


double antithetic_variates(double n)
		{
			//seed rand function
		 srand48(time(NULL)); //Make a new seed for the random number generator
		 double x=0, y=0, S_n=0, errorterm=0, error=0, g=0, y_comp=0;
			// loop according to sample size n
			for (int i = 0; i < n; i++) {

					x = drand48();
					//x = 2-sqrt(4-3*z);
					y = (4/(1+x*x));
					y_comp = (4/(1+(1-x)*(1-x)));

					S_n += (y + y_comp)/2;
					errorterm += (((y + y_comp)/2)*((y + y_comp)/2));
			}
			S_n = S_n/n;
			error = sqrt(errorterm/n - S_n*S_n)/sqrt(n);
		printf("Antithetic Variates:\n");
		printf("S_n= %f, Estimated error +/- %f, Real error =%f\n", S_n, error, (M_PI-S_n));
		}


      double stratified_sampling(double n)
					{
						//seed rand function
					 srand48(time(NULL)); //Make a new seed for the random number generator
					int k = 4; // Number of domains
					double x=0, y=0, S_n=0, errorterm=0, error=0;
						// loop according to sample size n
			    		double M_j = (double)1/(double)4; // domain size
				            int nj = n/k; // domain sample size
				                // loop according to number of domains
				                for (int i = 1; i < k+1; i++)
				                {
				                    double S_n1=0, S_n2=0, sigma=0;
				                        for (int j = 0; j < nj; j++)
				                        {

								            x = drand48()/(double)k + (double)(i-1)/(double)k; // uniformly distributed
								            y = (4/(1+x*x));
								            S_n1 += y*y;
								            S_n2 += y;

				                        }
								   S_n += M_j*(S_n2/nj);
								   S_n1 = S_n1/nj;
								   S_n2 = S_n2/nj;
								   sigma = (S_n1 - (S_n2*S_n2));
								   errorterm += (M_j*M_j)*(sigma/nj);
						        }

						error = sqrt(errorterm);
					printf("Stratified Sampling:\n");
					printf("S_n= %f, Estimated error +/- %f, Real error =%f\n", S_n, error, (M_PI-S_n));
					}


int main() {

	    for (int exp = 5; exp <8; exp++)
	    {
	    double size = pow(10,exp);
	    printf("For n = 10^%d-------------------------------------\n", exp);

		ordinary_montecarlo(size);
	    importance_sampling(size);
	    control_variates(size);
	    antithetic_variates(size);
	    stratified_sampling(size);
        }

return 0;
}
