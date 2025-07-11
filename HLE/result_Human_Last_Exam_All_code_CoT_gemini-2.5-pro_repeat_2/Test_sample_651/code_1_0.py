import math

# Based on the analytical derivation, the limit of M(theta) is calculated.
# The expression for M(theta) approaches the following form as theta -> 0:
# lim M(theta) = lim [pi/2 - atan((2/3)*cot(theta/2))]

# As theta approaches 0 from the positive side:
# 1. cot(theta/2) approaches positive infinity.
# 2. The argument of atan, (2/3)*cot(theta/2), also approaches positive infinity.
# 3. The arctan of a value approaching positive infinity is pi/2.

# This leads to the final equation for the limit calculation.
pi_value = math.pi
limit_of_atan_term = pi_value / 2
limit_of_M = (pi_value / 2) - limit_of_atan_term

# The problem asks to output each number in the final equation.
# The final equation is: (pi/2) - (pi/2) = 0.
# The numbers are shown below.
print("The final calculation for the limit is based on the equation:")
print("lim M(theta) = pi/2 - lim atan(...)")
print("\nThe numbers in the final equation are:")
print(f"First term (pi/2): {pi_value / 2}")
print(f"Second term (limit of atan): {limit_of_atan_term}")
print(f"Result of the subtraction: {limit_of_M}")
print(f"\nThus, the final equation is: {pi_value / 2} - {limit_of_atan_term} = {limit_of_M}")
