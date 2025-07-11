import math

# The goal is to find the minimum q such that Alice's win probability p_A is >= 0.95.
# This is equivalent to finding q when p_A is exactly 0.95.
p_A = 0.95

# The final equation we derived is q = [ (1 - (1 - p_A)^(1/3)) / p_A^3 ]^(1/4)
# We will calculate the result step by step, showing the numbers in the equation.

# 1. Calculate the term (1 - p_A)
val1 = 1 - p_A
print(f"The value of (1 - p_A) is: {val1}")

# 2. Calculate its cube root: (1 - p_A)^(1/3)
val2 = val1**(1/3)
print(f"The value of (1 - p_A)^(1/3) is: {val2}")

# 3. Calculate the numerator: 1 - (1 - p_A)^(1/3)
numerator = 1 - val2
print(f"The numerator 1 - (1 - p_A)^(1/3) is: {numerator}")

# 4. Calculate the denominator: p_A^3
denominator = p_A**3
print(f"The denominator p_A^3 is: {denominator}")

# 5. Calculate q^4
q_fourth_power = numerator / denominator
print(f"The value of q^4 is: {q_fourth_power}")

# 6. Calculate q (q_0 in the problem description)
q0 = q_fourth_power**(1/4)
print(f"The minimum value of q (q0) is: {q0}")

# 7. Calculate floor(100 * q0)
result = math.floor(100 * q0)
print(f"\nThe final result for floor(100 * q0) is:")
print(result)