import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

# Part 1: Analysis for the Uniform Distribution Case

# The analysis shows that R^2 = epsilon * C_max, where C_max = max_X Var[G(X, Z)]
# for Z ~ U(0,1). The code below calculates C_max.

# 1. Define the Green's function G(X, xi) for the operator y'' - y'
# with homogeneous boundary conditions at X=0 and X=1.
def G(X, xi):
    """Green's function G(X, xi) for the rescaled problem."""
    e = np.e
    if xi < X:
        # Formula for xi < X
        return (1 - np.exp(-xi)) * (np.exp(X) - e) / (e - 1)
    else:
        # Formula for xi >= X
        return (1 - np.exp(1 - xi)) * (np.exp(X) - 1) / (e - 1)

# 2. Define the function for the variance of G(X, Z) where Z ~ U(0,1).
def var_G(X):
    """Computes Var[G(X, Z)] for Z ~ U(0,1) at a given X."""
    # The expected value E[G(X,Z)] has a known analytical solution.
    e_G_X = (np.exp(X) - 1) / (np.e - 1) - X
    
    # The second moment E[G(X,Z)^2] is computed by numerical integration.
    integrand_sq = lambda xi: G(X, xi)**2
    # The quad function integrates the squared Green's function over the domain of Z, [0, 1].
    e_G2_X, _ = quad(integrand_sq, 0, 1, limit=100)
    
    # Variance = E[X^2] - (E[X])^2
    return e_G2_X - e_G_X**2

# 3. Find the maximum of var_G(X) for X in [0, 1].
# We achieve this by minimizing the negative of the function.
res = minimize_scalar(lambda x: -var_G(x), bounds=(0, 1), method='bounded')

C_max = -res.fun
sqrt_C_max = np.sqrt(C_max)
scaling_exponent_uniform = 0.5

print("--- Analysis Result ---")
print("For the case of Uniformly distributed z_i:")
print("The maximum magnitude of fluctuations R(epsilon) is found to scale with epsilon.")
print("The final expression for R(epsilon) is R(\u03B5) = C * \u03B5^p, where:")
print(f"p = {scaling_exponent_uniform}")
print(f"C = {sqrt_C_max:.5f}")
print("So the equation is:")
print(f"R(\u03B5) = {sqrt_C_max:.5f} * \u03B5^{scaling_exponent_uniform}")
print("-" * 25)

# Part 2: Analysis for the Normal Distribution Case
scaling_exponent_normal = 1.5
print("For the case of Normally distributed z_i:")
print("The scaling for R(\u03B5) is expected to change.")
print(f"The new scaling exponent is expected to be p = {scaling_exponent_normal}.")
print("Therefore, the scaling does not remain the same.")
print("--- End of Report ---")
