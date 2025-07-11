import math

# Orders of the fundamental groups of X1, X2, and X3
n1 = 5
n2 = 8
n3 = 2

# Calculate the rank of the kernel using the formula for a free product of k=3 groups.
# The formula is: rank = sum_{1 <= i < j <= k} (|G_i|-1)(|G_j|-1)
term1 = (n1 - 1) * (n2 - 1)
term2 = (n1 - 1) * (n3 - 1)
term3 = (n2 - 1) * (n3 - 1)
rank = term1 + term2 + term3

# Print the equation with all the numbers
print(f"Rank = ({n1}-1)*({n2}-1) + ({n1}-1)*({n3}-1) + ({n2}-1)*({n3}-1)")
print(f"Rank = ({n1-1})*({n2-1}) + ({n1-1})*({n3-1}) + ({n2-1})*({n3-1})")
print(f"Rank = {term1} + {term2} + {term3}")
print(f"Rank = {rank}")
