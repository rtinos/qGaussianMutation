# qGaussianMutation
The code is used to produce q-Gaussian mutations. The shape of the q-Gaussian mutation distribution is controlled by a real parameter q.

************** q-Gaussian Mutation **************

Reference: Tinos & Yang (2011). "Use of the q-Gaussian mutation in evolutionary algorithms", Soft Computing, 15(8): 1523-1549.

The code has 2 functions that create an offspring from a parent (represented by a vector with lcrom real numbers). The functions are:

1) void mutation_q_gaussian_isotropic( double *parent , double *offspring , int lcrom, double *mut_parameter, double q ) 
uses an isotropic q-Gaussian distribution with parameter q. 

2) void mutation_q_gaussian_anisotropic( double *parent , double *offspring , int lcrom, double *mut_parameter, double q ) 
uses an anisotropic q-Gaussian distribution with parameter q. 

Observation: 
a) the parameter mutation vector mut_parameter is similar to the vector of mutation parameters sigma in mutation generated from a Gaussian distribution. 
b) function gasdev generates a standard Gaussian random variable (scalar) - see book Numerical recipes in C
