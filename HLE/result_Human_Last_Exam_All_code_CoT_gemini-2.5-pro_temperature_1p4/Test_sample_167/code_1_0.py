# Set the parameter n for the graph K_n,n
n = 1000

# For a complete bipartite graph K_n,n where n is even,
# the Alon-Tarsi number is given by the formula: AT(K_n,n) = n/2 + 2.

# Perform the calculation based on the formula.
# First, divide n by 2.
half_n = n / 2
# Then, add 2 to the result.
result = half_n + 2

# Print the explanation of the calculation.
print(f"The Alon-Tarsi number of K_{n},{n} for an even n is calculated using the formula n/2 + 2.")
print("The calculation steps are:")
print(f"1. Divide n by 2: {n} / 2 = {int(half_n)}")
print(f"2. Add 2 to the result: {int(half_n)} + 2 = {int(result)}")
print(f"\nFinal Equation: {n} / 2 + 2 = {int(result)}")
print(f"The Alon-Tarsi number of K_{n},{n} is {int(result)}.")