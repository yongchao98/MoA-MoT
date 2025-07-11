# Let C be a genus 3 Riemann surface and X its 15th symmetric power.
# This script calculates the smallest and the largest possible rank of the Neron-Severi group of X.

# Define the parameters based on the problem statement
g = 3  # genus of the curve C
d = 15 # degree of the symmetric power X

# The rank of the Neron-Severi group of X, rho(X), is related to the rank
# of the Neron-Severi group of the Jacobian of C, rho(Jac(C)).
# The condition d >= 2g - 1 is met (15 >= 5), so we can use the formula:
# rho(X) = rho(Jac(C)) + 1

# --- Smallest Rank Calculation ---
# The smallest possible rank for rho(Jac(C)) of a genus g=3 curve is 1.
min_rho_jac = 1
smallest_rank_X = min_rho_jac + 1

print("Smallest possible rank of the Neron-Severi group of X:")
print("rho_min(X) = rho_min(Jac(C)) + 1 = {} + 1 = {}".format(min_rho_jac, smallest_rank_X))
print("-" * 20)

# --- Largest Rank Calculation ---
# The largest possible rank for rho(Jac(C)) of a genus g=3 curve is g^2.
max_rho_jac = g**2
largest_rank_X = max_rho_jac + 1

print("Largest possible rank of the Neron-Severi group of X:")
print("rho_max(X) = rho_max(Jac(C)) + 1 = {} + 1 = {}".format(max_rho_jac, largest_rank_X))

print("\nTherefore, the smallest possible rank is 2 and the largest is 10.")