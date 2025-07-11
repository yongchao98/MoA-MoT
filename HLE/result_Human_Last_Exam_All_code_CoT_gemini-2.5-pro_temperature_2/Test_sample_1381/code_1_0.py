# The number of species in the system.
# This value can be changed to calculate the number of equilibria for a different system size.
N = 10

# The maximum number of equilibria is 2^N, which can be calculated using the exponentiation operator.
base = 2
exponent = N
num_equilibria = base ** exponent

# Print the result in a clear format that shows the calculation,
# satisfying the requirement to output each number in the final equation.
print(f"For a system with N = {N} species, the maximum possible number of equilibria is:")
print(f"{base}^{exponent} = {num_equilibria}")