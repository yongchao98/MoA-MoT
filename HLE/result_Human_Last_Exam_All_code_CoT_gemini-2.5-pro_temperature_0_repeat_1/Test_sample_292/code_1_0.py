import math

# The problem is to calculate the sum S for n=99.
n = 99

# Based on mathematical derivation, the sum S simplifies to the formula:
# S = (n + 1)^(n - 1)

# We substitute n = 99 into the formula.
base = n + 1
exponent_of_base = n - 1

# The result is 100^98. We need to express this as a power of 10.
# 100 = 10^2
# So, S = (10^2)^98 = 10^(2 * 98)
final_base = 10
final_exponent = 2 * exponent_of_base

# As requested, we output each number in the final equation.
print(f"The simplified formula for the sum is S = (n + 1)^(n - 1).")
print(f"For n = {n}, the equation is:")
print(f"S = ({n} + 1)^({n} - 1)")
print(f"S = {base}^{exponent_of_base}")
print(f"Expressing this as a power of 10:")
print(f"S = ({final_base}^2)^{exponent_of_base}")
print(f"S = {final_base}^(2 * {exponent_of_base})")
print(f"S = {final_base}^{final_exponent}")