import numpy as np
from scipy.integrate import quad

# Step 1: Define the constants
phi = (1 + np.sqrt(5)) / 2
e = np.e

# Step 2: Calculate the upper integration limit
upper_limit = phi**3 - 1

# Step 3: Implement the simplified integrand function
def integrand(x):
    """
    Calculates the value of the integrand:
    cos(ln(1 + exp(arctan(ln(cos(x/e))))))
    """
    # The argument to log is cos(x/e). For the given limits, this is > 0.
    # To prevent potential floating point issues near the boundaries of a function's domain,
    # we can add a small epsilon or clip the value, but for the given range,
    # cos(x/e) is safely away from zero.
    cos_val = np.cos(x / e)
    
    # Chain the functions as in the expression
    log_cos_val = np.log(cos_val)
    arctan_val = np.arctan(log_cos_val)
    exp_val = np.exp(arctan_val)
    final_arg = np.log(1 + exp_val)
    
    return np.cos(final_arg)

# Step 4: Use scipy.integrate.quad to compute the definite integral
# quad returns a tuple (result, estimated_error)
integral_value, error = quad(integrand, 0, upper_limit)

# Step 5: Print the final result. The problem statement asks to
# "output each number in the final equation". The most important number
# in the equation `Integral = Value` is the final value.
print(integral_value)
