# Values of the exponents derived from the physical analysis of the system.
n1 = 1
n2 = 1
n3 = -1
n4 = 0
n5 = -2
n6 = -1.5

# The sum to be calculated is sum_{k=1 to 6} k * n_k
# This corresponds to: 1*n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6
total_sum = 1*n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6

# We print the equation with the numerical values
# The formatting is to ensure that signs are handled correctly for display
print(f"The equation is: 1*({n1}) + 2*({n2}) + 3*({n3}) + 4*({n4}) + 5*({n5}) + 6*({n6})")

# Calculate each term
term1 = 1 * n1
term2 = 2 * n2
term3 = 3 * n3
term4 = 4 * n4
term5 = 5 * n5
term6 = 6 * n6

# We also print the evaluation of the equation
print(f"Evaluation: {term1} + {term2} + {term3} + {term4} + {term5} + {term6} = {total_sum}")

# Finally, we print the result of the calculation
print("\nThe value of the sum is:")
print(total_sum)