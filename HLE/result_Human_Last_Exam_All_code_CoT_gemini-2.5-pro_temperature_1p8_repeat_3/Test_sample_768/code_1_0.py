import numpy as np
from scipy.integrate import quad

# 1. Define the constants and integration limits
phi = (1 + np.sqrt(5)) / 2
lower_limit = 0
upper_limit = phi**3 - 1

# 2. Define the integrand function based on the simplification
def integrand(x):
    """
    Calculates the value of the integrand:
    Re[1 / (1 + exp(arctan(log(cos(x/e)))))**i]
    which simplifies to cos(log(1 + exp(arctan(log(cos(x/e))))))
    """
    # To prevent potential floating-point issues at x=0, we can handle it as a special case,
    # though numpy handles log(1.0) correctly.
    # At x=0, cos(0)=1, log(1)=0, arctan(0)=0, exp(0)=1, result is cos(log(2))
    if x == 0:
        return np.cos(np.log(2.0))
    
    val_cos = np.cos(x / np.e)
    # The arguments to log, arctan, exp are all well-defined in the integration range
    val_log_inner = np.log(val_cos)
    val_atan = np.arctan(val_log_inner)
    val_exp = np.exp(val_atan)
    
    return np.cos(np.log(1 + val_exp))

# 3. Perform the numerical integration
result, error = quad(integrand, lower_limit, upper_limit)

# 4. Output the numbers related to the evaluation of the integral.
# The integral is of the form: integral from A to B of f(x) dx = C
# The numbers in this equation are A, B, and C.
print(f"The final equation is: integral from {lower_limit} to {upper_limit} = {result}")
print("The numbers in this equation are:")
print(f"Lower limit: {lower_limit}")
print(f"Upper limit: {upper_limit}")
print(f"Value of the integral: {result}")
