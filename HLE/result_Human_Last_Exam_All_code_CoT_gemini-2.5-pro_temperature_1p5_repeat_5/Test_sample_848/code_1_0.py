import math

# The problem is to find the integer part of 10^4 * L, where L is the limit of F(N)/ln(N).
# Based on the analysis, the solutions (a,b) fall into two families,
# characterized by an invariant k=13 or k=5.
# The growth rate of the solution sequence in each family is the largest root of r^2 - k*r + 1 = 0.

# For k=13, the growth rate is r_A.
k_A = 13
r_A = (k_A + math.sqrt(k_A**2 - 4)) / 2

# For k=5, the growth rate is r_B.
k_B = 5
r_B = (k_B + math.sqrt(k_B**2 - 4)) / 2

# The limit L is given by the sum of the contributions from each family.
# The factor of 2 accounts for symmetric pairs (a,b) and (b,a).
L = 2 / math.log(r_A) + 2 / math.log(r_B)

# The final quantity to compute.
final_value = 10000 * L

# Output the numbers used in the final equation.
print(f"The equation for the limit L is: 2/ln(r_A) + 2/ln(r_B)")
print(f"Value for r_A (k={k_A}): {r_A}")
print(f"Value for r_B (k={k_B}): {r_B}")
print(f"The limit L is: {L}")
print(f"The expression to be computed is: 10000 * L = {final_value}")

# The problem asks for the integer part of the final value.
integer_part = math.floor(final_value)
print(f"The integer part of 10000 * L is: {integer_part}")

# Return the answer in the specified format
# The final result is the integer part of 10^4 * L
# answer = math.floor(10000 * (2/math.log((13+math.sqrt(165))/2) + 2/math.log((5+math.sqrt(21))/2)))
# <<<20580>>>