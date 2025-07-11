#
# Description:
# This script calculates the trace of the covariance matrix for the random variable v
# based on the analytical derivation outlined in the problem description.
# The final formula is Tr(Cov(v)) = 1 - ((alpha - beta) / (alpha + beta))^2.
#

# Given parameters
alpha = 3
beta = 2

# Step 1: Calculate the numerator and denominator of E[d_1]
# E[d_1] = (alpha - beta) / (alpha + beta)
numerator_E_d1 = alpha - beta
denominator_E_d1 = alpha + beta

# Step 2: Calculate E[d_1]
E_d1 = numerator_E_d1 / denominator_E_d1

# Step 3: Calculate the squared norm of the expectation of d, which is (E[d_1])^2
norm_E_d_squared = E_d1 ** 2

# Step 4: Calculate the trace of the covariance matrix
# Tr(Cov(v)) = 1 - ||E[d]||^2
trace_covariance = 1 - norm_E_d_squared

# Output the explanation and the final equation with all numbers.
print("The calculation of the trace of the covariance matrix follows the formula:")
print("Trace = 1 - ( (alpha - beta) / (alpha + beta) )^2")
print("\nSubstituting the given values alpha = {} and beta = {}:".format(alpha, beta))
print("Trace = 1 - ( ({} - {}) / ({} + {}) )^2".format(alpha, beta, alpha, beta))
print("Trace = 1 - ( {} / {} )^2".format(numerator_E_d1, denominator_E_d1))
print("Trace = 1 - ( {} )^2".format(E_d1))
print("Trace = 1 - {}".format(norm_E_d_squared))
print("Trace = {}".format(trace_covariance))