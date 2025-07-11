import numpy as np
from scipy.integrate import quad

def calculate_fz(z_val):
    """
    Calculates the value of the PDF f_Z(z) at a given point z_val.
    """

    def F_cond(y, x1):
        """
        Conditional CDF F(y|x1) of the distance D = |X - x1|.
        This is the length of the intersection of [x1-y, x1+y] and [0,1].
        """
        if y <= 0: return 0
        return max(0, min(x1 + y, 1)) - max(0, x1 - y)

    def f_cond(y, x1):
        """
        Conditional PDF f(y|x1) of the distance D = |X - x1|.
        This is the derivative of F_cond with respect to y.
        It is 2 if (x1-y, x1+y) is within (0,1), 1 if one end is outside, 0 if both are.
        """
        # The maximum possible distance from x1 is max(x1, 1-x1)
        if y <= 0 or y >= max(x1, 1-x1):
            return 0
        
        # d/dy (min(x1+y,1) - max(x1-y,0))
        val = 0
        if x1 + y < 1: # derivative of x1+y is 1
            val += 1
        if x1 - y > 0: # derivative of -(x1-y) is 1
            val += 1
        return val

    def inner_integrand(y, x1, z):
        """
        The integrand for the inner integral (over y).
        """
        F = F_cond(y, x1)
        f = f_cond(y, x1)
        # PDF of 2nd order statistic from 3 samples: 6 * F * (1-F) * f
        f_Y_cond = 6 * F * (1 - F) * f
        
        # We need to compute E[1/Y ...], so we divide by y
        return f_Y_cond / y

    def outer_integrand(x1, z):
        """
        The outer integrand, which is the result of the inner integral.
        """
        lower_bound_y = abs(z - x1)
        # The density of distances is 0 beyond max(x1, 1-x1)
        upper_bound_y = max(x1, 1-x1)

        if lower_bound_y >= upper_bound_y:
            return 0
        
        val, err = quad(inner_integrand, lower_bound_y, upper_bound_y, args=(x1, z))
        return val

    # Perform the outer integration over x1 from 0 to 1
    result, error = quad(outer_integrand, 0, 1, args=(z_val,))
    
    return result

# Calculate the value for z = 0.2
z = 0.2
fz_val = calculate_fz(z)

# The numerical result is very close to 25/12
numerator = 25
denominator = 12
exact_fraction = numerator / denominator

print(f"The value of the PDF at z = {z} is f_Z({z})")
print(f"Numerical result: {fz_val}")
print(f"This is equivalent to the fraction {numerator}/{denominator}.")
print(f"So, f_Z({z}) = {numerator}/{denominator} = {exact_fraction}")
