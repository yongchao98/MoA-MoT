import sys
from fractions import Fraction

# This script calculates the value of c based on a key result from harmonic analysis.

# Step 1: Define the given Frostman dimension alpha.
# The problem states we have a 8/5-Frostman measure.
alpha_num = 8
alpha_den = 5
alpha = alpha_num / alpha_den

# Step 2: State the formula for the exponent c.
# The problem asks for the smallest possible c such that the L2 norm of hat(mu)(r*sigma)
# is O(r^(c+epsilon)). This corresponds to finding the sharp decay exponent.
# A result by T. Wolff gives that for a measure on R^2 with Frostman dimension alpha > 3/2,
# the integral I(r) = integral_{S^1} |hat(mu)(r*sigma)|^2 d(sigma) is O(r^p)
# with p = -(2*alpha - 3)/2.
# The L2 norm is sqrt(I(r)), so its exponent is p/2.
# Therefore, c = p/2 = -(2*alpha - 3)/4.

# We must check if alpha > 3/2.
# alpha = 8/5 = 1.6. 3/2 = 1.5. Since 1.6 > 1.5, the formula is applicable.

# Step 3: Calculate c using the formula.
# We will show the calculation step-by-step.
print("The formula for the exponent c is: c = -(2 * alpha - 3) / 4")
print(f"Given alpha = {alpha_num}/{alpha_den}")
print("\nSubstituting the value of alpha into the formula:")
print(f"c = -(2 * ({alpha_num}/{alpha_den}) - 3) / 4")
print(f"c = -({2 * alpha_num}/{alpha_den} - 3) / 4")
# To subtract, we find a common denominator
print(f"c = -({2 * alpha_num}/{alpha_den} - {3 * alpha_den}/{alpha_den}) / 4")
numerator = 2 * alpha_num - 3 * alpha_den
print(f"c = -({numerator}/{alpha_den}) / 4")
print(f"c = -{numerator} / ({alpha_den} * 4)")

final_num = -numerator
final_den = alpha_den * 4
c = final_num / final_den

# Step 4: Display the final result.
# The result can be expressed as a fraction and a decimal.
final_fraction = Fraction(c).limit_denominator()
print("\nThe final calculated value is:")
print(f"c = {final_num}/{final_den} = {final_fraction}")
print(f"As a decimal, c = {c}")

# To conform to the output format, we also output the final answer separately.
# Python 2/3 compatibility for the final answer format
if sys.version_info[0] < 3:
    # In Python 2, we might not get a float from division
    c_final_float = float(final_num)/final_den
    print("<<<{}>>>".format(c_final_float))
else:
    print("<<<{}>>>".format(c))