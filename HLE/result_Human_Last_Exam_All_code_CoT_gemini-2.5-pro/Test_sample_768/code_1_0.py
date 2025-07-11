import numpy as np
from scipy.integrate import quad

# Step 1: Define the constants and limits of integration
# Golden ratio
phi = (1 + np.sqrt(5)) / 2

# Lower limit of integration
lower_limit = 0

# Upper limit of integration is phi^3 - 1, which simplifies to 2*phi
upper_limit = 2 * phi

# Step 2: Define the integrand function
# The integrand is Re[1 / (1 + exp(arctan(log(cos(x/e)))))**i]
# which simplifies to cos(log(1 + exp(arctan(log(cos(x/e))))))
def integrand(x):
    # To avoid errors, ensure the argument of log is positive.
    # For the given range, cos(x/np.e) is always positive.
    cos_val = np.cos(x / np.e)
    if cos_val <= 0:
        return 0 # Should not happen in the integration interval
    
    inner_log = np.log(cos_val)
    arctan_val = np.arctan(inner_log)
    exp_val = np.exp(arctan_val)
    outer_log = np.log(1 + exp_val)
    return np.cos(outer_log)

# Step 3: Perform the numerical integration
# quad returns the result and an estimate of the absolute error
result, error = quad(integrand, lower_limit, upper_limit)

# Step 4: Output the result as requested by the final equation: Integral = value
print(f"The value of the integral from {lower_limit} to {upper_limit:.4f} is:")
print(result)
