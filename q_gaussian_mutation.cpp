/*******************************************************************************************************************************\
*				 q-Gaussian Mutation			   							*
*									   							*
* See Tinos & Yang (2011). "Use of the q-Gaussian mutation in evolutionary algorithms", Soft Computing, 15(8): 1523-1549	*
*																*
* mutation_q_gaussian_isotropic(): Isotropic q-Gaussian mutation 								*
* mutation_q_gaussian_anisotropic(): Anisotropic q-Gaussian mutation 								*
*				Inputs: parent solution, parameter q, size of the solutions (lcrom), mutation parameter vector 	*
*				Output: offspring solution									*
*																*
\*******************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/***********************************************************************\
*	q-Exponential of x						*
\***********************************************************************/
double log_q( double x, double q ){

	if (fabs(q - 1.0) < EPS)
		return(log(x));					// If q is 1, use the usual natural logarithm
	else
		return( (pow(x,1.0-q) - 1.0)/(1.0-q));		// If q differs from 1, use the definition of the q-log

}


/*******************************************************************************\
*   Generate a random variable from a q-Gaussian Distribution (Tsallis)		*
* see Thistleton, Marsh, Nelson & Tsallis (2006)				*							
*							                        *
* Returns random deviates drawn from a q-gaussian distribution			*
* The q that characterizes the q-Gaussian is given by qDist		        *
* Obs.: q must be < 3							        *
\*******************************************************************************/
double dist_q_gaussian( double qDist ){

	double qGen, u1, u2;

	// Compute the q to be used on the q-log
	qGen = (1.0 + qDist)/(3.0 - qDist);

	// Get two uniform random deviates
 	u1=rand()/(RAND_MAX+1.0);
	u2=rand()/(RAND_MAX+1.0);
	
	// Apply the generalized Box-Muller algorithm
	return( sqrt(2*fabs(log_q(u1,qGen)))*sin(2*PI*u2) );

}


/***********************************************************************************\
*	Euclidean norm of x	  	    				    	    *
\***********************************************************************************/
double norm_euc(double *x, int l)
{
	double norm=0.0;
	int i;

	for (i=0;i<l;i++)
		norm += x[i]*x[i];

	return ( sqrt(norm) );
}


/***********************************************************************************\
*	Isotropic q-Gaussian mutation   	    				    *
\***********************************************************************************/
void mutation_q_gaussian_isotropic( double *parent , double *offspring , int lcrom, double *mut_parameter, double q ){

	int i;
	double norm_vector, norm_vector_new;

	norm_vector_new = dist_q_gaussian(q);	// norm of the vector given by a random variable from q-Gaussian distribution
	for (i=0;i<lcrom;i++)
	  offspring[i] = gasdev(&idum); 	// random variable from Gaussian distribution 
					   	// Here, the function gasdev was used 
						// (see book Numerical recipes in C)					      
	norm_vector = norm_euc(offspring,lcrom);	// returns the Euclidean norm of the vector
	for (i=0;i<lcrom;i++){
		offspring[i] = parent[i] + mut_parameter[i]*norm_vector_new*offspring[i]/norm_vector ;
		// Obs.: if necessary, check the limits of offspring[i]		
	}

}


/***********************************************************************************\
*	Anisotropic q-Gaussian mutation					  	    *
\***********************************************************************************/
void mutation_q_gaussian_anisotropic( double *parent , double *offspring , int lcrom, double *mut_parameter, double q  ){

	int i;

	for (i=0;i<lcrom;i++){
		offspring[i] = parent[i] + mut_parameter[i]*dist_q_gaussian(q);
		// Obs.: if necessary, check the limits of offspring[i]		
	}

}
