import math

# The user can change this value for n.
n = 10

# The approximation formula is S_n ≈ n^2/2 + 1/120 + 1/(252*n)
term1 = n**2 / 2
term2 = 1 / 120
term3 = 1 / (252 * n)

# Calculate the sum of the terms
result = term1 + term2 + term3

# The problem requires outputting each number in the final equation.
# The following print statement displays the formula with the specific value of n and the numbers.
print(f"The approximation for n = {n} is derived from the formula: n^2/2 + 1/120 + 1/(252*n)")
print(f"Calculation: ({n}^2)/2 + 1/120 + 1/(252*{n}) = {term1} + {term2:.6f} + {term3:.6f}")
print(f"The final approximated value is: {result}")
