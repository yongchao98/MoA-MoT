# Define the properties of the surface Sigma
g_sigma = 10  # Genus of Sigma
b_sigma = 1   # Number of boundary components of Sigma

# The smallest positive integer g for the containing closed surface Sigma'
# is determined by the genus of the "double" of Sigma. This represents
# the upper bound that works for all possible embeddings, and it can be
# shown to be the lower bound for a "worst-case" embedding.
# The formula is: g = 2 * g_sigma + b_sigma - 1.

# Calculate the genus g
g = 2 * g_sigma + b_sigma - 1

# Output the final equation with the numbers plugged in
print(f"Given a surface Sigma with genus = {g_sigma} and number of boundary components = {b_sigma}.")
print("The smallest genus g of a closed surface Sigma' that must contain Sigma is given by the formula:")
print("g = 2 * (genus of Sigma) + (number of boundary components) - 1")
print(f"g = 2 * {g_sigma} + {b_sigma} - 1")
print(f"g = {g}")