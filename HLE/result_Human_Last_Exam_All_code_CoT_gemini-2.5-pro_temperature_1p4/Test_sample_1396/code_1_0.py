import math

# The problem asks for a realistic upper bound for a connected epsilon-envy-free
# allocation for n=4 agents. We refer to the paper cited in the prompt,
# "The query complexity of connected-piece cake cutting" by Brânzei and Nisan (2022).
# Their algorithm for n agents has a query complexity of O(n! / epsilon^(n-1)).
# We will calculate the constant factor n! for n=4.

n = 4

# Calculate factorial of n
n_factorial = math.factorial(n)

# Print the final equation as requested
# 4! = 1 * 2 * 3 * 4 = 24
numbers = list(range(1, n + 1))
equation_str = " * ".join(map(str, numbers))
print(f"For n = {n} agents, the constant factor in the upper bound from the cited paper is based on the factorial of n.")
print(f"The equation is: n! = {equation_str}")
print(f"The result is: {n_factorial}")
