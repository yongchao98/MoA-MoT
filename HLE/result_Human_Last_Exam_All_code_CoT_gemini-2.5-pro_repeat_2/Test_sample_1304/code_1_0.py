import numpy as np
from scipy.optimize import minimize_scalar

# We found that the problem reduces to maximizing I(z1, 1) = -5*z1*(z1+1) for z1 in [-1, 0].
# To use a minimizer, we find the minimum of the negative of this function.
def neg_I_z2_is_1(z1):
    """The negative of the function I(z1, 1) to be minimized."""
    return -(-5 * z1 * (z1 + 1))

# We use scipy.optimize.minimize_scalar to find the z1 that maximizes I(z1, 1).
# The search is bounded to the interval [-1, 0].
result = minimize_scalar(neg_I_z2_is_1, bounds=(-1, 0), method='bounded')

# The optimal z1 and the maximum value of the integral I
optimal_z1 = result.x
max_I = -result.fun

# The parameters for the final calculation of c_3
c3_factor_num = 7
c3_factor_den = 2
max_I_num = 5
max_I_den = 4

# Calculate the maximum value of c_3
max_c3 = (c3_factor_num / c3_factor_den) * max_I

# The numerator and denominator of the final result
final_num = c3_factor_num * max_I_num
final_den = c3_factor_den * max_I_den

print("The problem of maximizing c_3 was reduced to maximizing an integral I.")
print(f"The maximum value of I was found by optimizing with respect to z1 and z2.")
print(f"The optimal value for z1 is: {optimal_z1:.4f}")
print(f"The corresponding maximum value of the integral I is: {max_I:.4f} which is {max_I_num}/{max_I_den}")
print("\nThe final equation for the maximum value of c_3 is:")
print(f"max c_3 = ({c3_factor_num}/{c3_factor_den}) * ({max_I_num}/{max_I_den}) = {final_num}/{final_den}")
print(f"\nThe maximum value of c_3 is {max_c3}")