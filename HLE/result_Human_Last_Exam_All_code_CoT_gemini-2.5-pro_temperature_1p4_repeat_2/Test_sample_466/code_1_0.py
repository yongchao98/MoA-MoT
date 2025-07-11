import math

q = 11
n = 3

# Calculate the terms
term1 = q**n - 1
term2 = q**n - q
term3 = q**n - q**2

# Calculate the total number of adjunctions
total_adjunctions = term1 * term2 * term3

# Print the final equation as requested
print(f"{term1} * {term2} * {term3} = {total_adjunctions}")