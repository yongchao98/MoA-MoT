# Based on the derivation, the upper bound for the expression can be formulated as
# ||B Q_{0, M}||_inf <= Factor * sqrt(N).
# Our goal is to find the maximum possible value for this Factor.

# The factor is given by 2 * (1 - beta_M).
# We determined that the value of beta_M is between 0 and 1 (inclusive).
# To find the upper bound for the factor, we need to find its maximum value.
# The expression 2 * (1 - beta_M) is maximized when beta_M is minimized.
min_beta_M = 0

# Calculate the maximum value of the factor.
constant_k = 2
max_factor = constant_k * (1 - min_beta_M)

# The final equation is: Factor = 2 * (1 - 0)
# We now print the numbers in this final equation, and the result.
print(f"The equation for the factor is: Factor = {constant_k} * (1 - {min_beta_M})")
print(f"The resulting upper-bound factor is: {max_factor}")

# Final Answer as per instructions
print("<<<2>>>")