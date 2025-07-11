import numpy as np
from scipy.special import iv
from scipy.optimize import brentq

# The problem is to find the largest x for which the summation equals 0.
# The summation is equivalent to the modified Bessel function of the first kind:
# S(x) = I_{x-1}(2)
# We need to solve the equation I_{x-1}(2) = 0.

# Let v = x - 1. We are looking for the largest root of the equation f(v) = I_v(2) = 0.
def bessel_func_of_order_v(v):
    """Represents the function I_v(2)"""
    return iv(v, 2)

# Based on the properties of Bessel functions, the largest root v_0 is known to be in the interval (-3, -2).
# Let's confirm there is a sign change in this interval.
# I_{-2}(2) = I_2(2) is positive.
# Let's check I_{-2.5}(2).
val_at_minus_2 = bessel_func_of_order_v(-2.0)
val_at_minus_2_5 = bessel_func_of_order_v(-2.5)

# We search for the root in the interval where the function changes sign.
# Since val_at_minus_2 > 0 and val_at_minus_2_5 < 0, the root is in (-2.5, -2.0).
lower_bound = -2.5
upper_bound = -2.0

# Use the Brent method to find the root v_0.
v_root = brentq(bessel_func_of_order_v, lower_bound, upper_bound)

# The problem asks for the value of x, and we have the relation v = x - 1.
# So, x = v + 1.
x_root = v_root + 1

# The final equation is I_{x-1}(2) = 0.
# Substituting the found x value gives I_{<x_root> - 1}(2) = 0, which is I_{<v_root>}(2) = 0.
# The question asks for the value of x.
print(f"{{{x_root:.3f}}}")