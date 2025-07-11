import numpy as np
from scipy.special import iv
from scipy.optimize import root_scalar

# The summation converges to the modified Bessel function of the first kind, I_{x-1}(2).
# We need to find the largest x for which I_{x-1}(2) = 0.
# Let v = x - 1. We are solving for the largest root of I_v(2) = 0.

# Define the function for which we want to find the root.
# The variable is 'v', and the parameter z is fixed at 2.
def bessel_func_for_root(v):
    return iv(v, 2.0)

# The roots of I_v(2) are all negative. We search for the largest root (closest to zero).
# By plotting or preliminary checks, we know the largest root is between -2.4 and -2.3.
# We use root_scalar with this bracket to find the root v.
sol = root_scalar(bessel_func_for_root, bracket=[-2.4, -2.3])
v_root = sol.root

# The value of x is v + 1.
x_val = v_root + 1

# As requested, output the numbers in the final equation I_v(z) = 0 and the resulting x.
# Here, z is 2. The equation solved is I_v(2) = 0.
print(f"The equation to solve is I_v(z) = 0 for v, with z=2.")
print(f"The largest root found for v is: {v_root:.5f}")
print(f"This corresponds to x = v + 1 = {x_val:.5f}")

# Print the final answer in the specified format {-a.bbb}
print("\nThe largest x value for which the summation converges to 0 is:")
print(f"{{{x_val:.3f}}}")