import numpy as np

# Given parameters
alpha = 3.0
beta = 2.0
d = 101.0

# The trace of the covariance matrix is given by the formula:
# Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2
# As derived in the explanation, E[||v||^2] = 1.
E_v_norm_sq_val = 1.0

# The squared norm of the expected value is ||E[v]||^2 = (E[d_1])^2
# E[d_1] = E[(a-b)/(a+b)] = 2 * E[a/(a+b)] - 1
# E[a/(a+b)] for a~Gamma(alpha,theta) and b~Gamma(beta,theta) is the mean of a Beta(alpha, beta) distribution.
# Mean of Beta(alpha, beta) = alpha / (alpha + beta)
mean_beta = alpha / (alpha + beta)

# Now calculate E[d_1]
E_d1 = 2 * mean_beta - 1

# And its square, which is ||E[v]||^2
norm_E_v_sq = E_d1**2

# Finally, calculate the trace
trace_covariance = E_v_norm_sq_val - norm_E_v_sq

# Print the step-by-step calculation with the given numbers
print("The final calculation for the trace of the covariance matrix is:")
print(f"Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2")
print(f"             = 1 - (2 * alpha/(alpha+beta) - 1)^2")
print(f"             = 1 - (2 * {int(alpha)}/({int(alpha)}+{int(beta)}) - 1)^2")
print(f"             = 1 - (2 * {mean_beta:.2f} - 1)^2")
print(f"             = 1 - ({E_d1:.2f})^2")
print(f"             = 1 - {norm_E_v_sq:.4f}")
print(f"             = {trace_covariance:.4f}")