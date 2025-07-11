import math

# The problem is to find 1/p_n for n=1000.
n = 1000

# The derived formula for 1/p_n is 4 * cos^2(pi / (n + 2)).
# The numbers in the final equation are 4, pi, and (n+2).
four = 4
denominator = n + 2

# Calculate the result.
result = four * (math.cos(math.pi / denominator))**2

# Output a clear description of the calculation and the final answer.
print(f"For n = {n}, the value of 1/p_n is given by the expression: 4 * cos^2(pi / (n + 2))")
print(f"Substituting n = {n}, we calculate: {four} * cos^2(pi / {denominator})")
print(f"The numerical value of 1/p_{n} is: {result}")