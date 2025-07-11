# Define parameters for the hypersurface X
n = 3  # Ambient space is CP^n
d = 5  # Degree of the hypersurface

# Step 1: Calculate the Euler characteristic of X
# Formula: chi = ((1-d)^(n+1) - 1)/d + (n+1)
chi = (pow(1 - d, n + 1) - 1) // d + (n + 1)
print(f"The Euler characteristic of X is chi = {chi}")

# Step 2: Calculate the second Betti number b2(X)
# For a smooth complex surface, chi = b0 - b1 + b2 - b3 + b4
# b0=1, b1=0, b3=0, b4=1. So, chi = 2 + b2.
b2 = chi - 2
print(f"The second Betti number of X is b2 = {b2}")

# Step 3: Determine the rank of the relative homotopy group pi_3(CP^3, X)
# This rank, k, is equal to b2(X) - 1.
k = b2 - 1
print(f"The rank of the relative homotopy group pi_3(CP^3, X) is k = b2 - 1 = {k}")

# Step 4: Calculate the rank of pi_3(X)
# The rank of pi_3(X) is equal to the rank of Gamma(Z^k), which is k*(k+1)/2.
rank_pi3_X = k * (k + 1) // 2
print(f"The rank of pi_3(X) is given by the formula k*(k+1)/2.")
print(f"rank(pi_3(X)) = {k} * ({k} + 1) / 2 = {rank_pi3_X}")
