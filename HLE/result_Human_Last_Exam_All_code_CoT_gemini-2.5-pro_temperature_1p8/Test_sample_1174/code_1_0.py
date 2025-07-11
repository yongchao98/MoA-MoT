# The exponents determined from the physical analysis
n1 = 1
n2 = 1
n3 = -2
n4 = 0
n5 = -2
n6 = -3/2

# The expression to be calculated is sum_{k=1 to 6} k * n_k
# Let's calculate each term of the sum
term1 = 1 * n1
term2 = 2 * n2
term3 = 3 * n3
term4 = 4 * n4
term5 = 5 * n5
term6 = 6 * n6

# Print the full equation for clarity
print("The calculation is:")
print(f"(1 * {n1}) + (2 * {n2}) + (3 * {n3}) + (4 * {n4}) + (5 * {n5}) + (6 * {n6})")
print(f"= {term1} + {term2} + {term3} + {term4} + {term5} + {term6}")

# Calculate the final sum
total_sum = term1 + term2 + term3 + term4 + term5 + term6

# Print the final result
print(f"The final sum is: {total_sum}")
