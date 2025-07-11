# The problem asks to calculate the sum of weights for all possible sequences.
# Through combinatorial analysis, the sum S_n for a vocabulary of size n
# and sequence length n is found to be S_n = (n+1)^(n-1).

# Set the value of n as given in the problem.
n = 99

# According to the derived formula, the base of the power is n + 1.
base = n + 1

# The exponent of the power is n - 1.
exponent = n - 1

# Calculate the final result using the simplified formula.
# Python's integers can handle arbitrarily large numbers, so we can compute this directly.
result = base ** exponent

# The problem asks to output the numbers in the final equation.
print(f"For n = {n}, the sum is calculated using the formula (n+1)^(n-1).")
print(f"The equation can be written as: ({n}+1)^({n}-1)")
print(f"Base of the power: {n} + 1 = {base}")
print(f"Exponent of the power: {n} - 1 = {exponent}")
print(f"This simplifies to: {base}^{exponent}")
print(f"As a power of 10, this is: (10^2)^{exponent} = 10^(2*{exponent}) = 10^{2*exponent}")
print(f"The final result is: {result}")