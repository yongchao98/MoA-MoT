import numpy as np
from scipy.special import iv
from scipy.optimize import fsolve

# This script finds the largest value of x for which the summation
# Sum_{i=0 to inf} 1 / ((x + i - 1)! * i!)
# converges to 0.

# As explained in the plan, the summation is equivalent to the modified Bessel
# function of the first kind, I_{x-1}(2). The problem is thus reduced
# to finding the largest root 'v' of I_v(2) = 0, where v = x - 1.

def bessel_iv_at_2(v):
    """
    This function represents I_v(2). We will find the root of this function.
    The input 'v' is the order of the Bessel function.
    """
    return iv(v, 2.0)

# The roots of I_v(2) for real v are all negative. The largest root (closest to 0)
# is known to be near -2.7. We use -2.5 as an initial guess for the numerical solver.
initial_guess = -2.5
v_root = fsolve(bessel_iv_at_2, initial_guess)[0]

# The value found is the largest root 'v'. We now calculate the corresponding 'x'
# using the relation x = v + 1.
x_value = v_root + 1

# The prompt requests that we output each number in the final equation.
# The final equation is x = v + 1. We display the values for v and x.
# The results are formatted to three decimal places.
print(f"The equation to solve is I_v(2) = 0, with v = x - 1.")
print(f"The largest root found is v = {v_root:.3f}")
print(f"The corresponding value of x is found from the equation: x = {v_root:.3f} + 1")
print(f"Therefore, the largest value for x is {x_value:.3f}")
print(f"\nThe final result in the requested format {-a.bbb} is:")
print(f"{{{x_value:.3f}}}")
