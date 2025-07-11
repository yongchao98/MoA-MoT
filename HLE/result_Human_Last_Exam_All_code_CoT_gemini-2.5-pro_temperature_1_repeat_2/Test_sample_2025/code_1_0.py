import numpy as np
from scipy import integrate

def F_D_given_x1(d, x1):
    """
    Calculates the CDF of the distance D = |X - x1|, where X ~ U[0,1].
    F(d) = P(D <= d | X1=x1)
    """
    if d < 0:
        return 0.0
    # The length of the interval [x1-d, x1+d] intersected with [0,1]
    return min(x1 + d, 1.0) - max(x1 - d, 0.0)

def f_D_given_x1(d, x1):
    """
    Calculates the PDF of the distance D = |X - x1|.
    f(d) is the derivative of F(d).
    It's 2 if d < min(x1, 1-x1), 1 if d is between min and max, and 0 otherwise.
    """
    if d < 0 or d > max(x1, 1 - x1):
        return 0.0
    val = 0
    if d < x1:
        val += 1
    if d < 1 - x1:
        val += 1
    return float(val)

def pdf_D2_given_x1(d, x1):
    """
    Calculates the PDF of the second order statistic D_(2) from a sample of 3 distances.
    f_D(2)(d) = 6 * F(d) * (1 - F(d)) * f(d)
    """
    if d < 0:
        return 0.0
    cdf_val = F_D_given_x1(d, x1)
    pdf_val = f_D_given_x1(d, x1)
    return 6.0 * cdf_val * (1.0 - cdf_val) * pdf_val

def inner_integrand(d, z, x1):
    """
    The integrand for the inner integral: (1/d) * f_{D_(2)|x1}(d|x1).
    """
    if d == 0:
        return 0.0 # This case will be handled by the integrator limits
    return (1.0 / d) * pdf_D2_given_x1(d, x1)

def f_Z_given_x1(x1, z):
    """
    Calculates the conditional PDF f_Z(z|x1) by integrating the inner_integrand.
    The integral is from |z-x1| to the maximum possible distance, max(x1, 1-x1).
    """
    d0 = abs(z - x1)
    upper_limit = max(x1, 1 - x1)
    if d0 >= upper_limit:
        return 0.0
    
    # We pass z and x1 as arguments to the integrand function
    val, err = integrate.quad(inner_integrand, d0, upper_limit, args=(z, x1))
    return val

def f_Z(z):
    """
    Calculates the final PDF f_Z(z) by integrating the conditional PDF f_Z(z|x1) over x1 from 0 to 1.
    """
    val, err = integrate.quad(f_Z_given_x1, 0, 1, args=(z,))
    return val

# The value we want to calculate
z_value = 0.2

# Perform the calculation
result = f_Z(z_value)

# Print the final equation and result
print(f"f(0.2) = {result}")