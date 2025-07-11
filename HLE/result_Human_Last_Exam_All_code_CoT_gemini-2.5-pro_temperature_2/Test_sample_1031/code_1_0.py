import math

# Define the grid size
n = 14

# The problem is equivalent to finding the minimum number of sites 
# needed to hit every 2x2 square in the grid.
# We use a constructive approach by placing infected sites in a diagonal pattern.
# A site (i, j) is infected if (i + j) mod 3 equals a chosen constant k.
# This pattern guarantees that every 2x2 square is hit.
# We need to find which choice of k (0, 1, or 2) yields the minimum number of sites.

# Initialize counters for each of the three possible residue classes modulo 3
counts = {0: 0, 1: 0, 2: 0}

# Iterate over every site in the n x n grid
for i in range(n):
    for j in range(n):
        # Increment the count for the corresponding residue class
        counts[(i + j) % 3] += 1

# Extract the counts for each class
sites_k0 = counts[0]
sites_k1 = counts[1]
sites_k2 = counts[2]

# The minimum number of sites for this method is the smallest of the three counts.
min_sites = min(sites_k0, sites_k1, sites_k2)

# Print the results of the calculation
print(f"To infect a {n}x{n} grid, we can infect all sites (i,j) where (i+j) mod 3 = k.")
print("The number of sites required for each choice of k is:")
print(f"k=0: Number of sites = {sites_k0}")
print(f"k=1: Number of sites = {sites_k1}")
print(f"k=2: Number of sites = {sites_k2}")
print(f"\nThe minimum number of sites needed using this construction is {min_sites}.")