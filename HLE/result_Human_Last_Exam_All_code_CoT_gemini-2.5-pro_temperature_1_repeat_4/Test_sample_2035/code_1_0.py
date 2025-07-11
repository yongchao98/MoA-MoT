import math

# Given parameters
n = 4048
# alpha_i^2 = i

# Based on the derivation, the expression simplifies to the sum of alpha_i^2 plus 10^15.
# We need to calculate the sum of the first n integers.
sum_alpha_sq = n * (n + 1) // 2

# The constant term to be added
constant_term = 10**15

# Calculate the final result
result = sum_alpha_sq + constant_term

# As requested, output the numbers in the final equation.
# The simplified final equation is: sum_alpha_sq + constant_term = result
print(f"The first term of the sum is: {sum_alpha_sq}")
print(f"The second term of the sum is: {constant_term}")
print(f"The final result of the sum is: {result}")