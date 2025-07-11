#
# This script calculates the smallest and largest possible ranks of the Neron-Severi
# group for the 15th symmetric power of a genus 3 Riemann surface.
#

# Genus of the Riemann surface C
g = 3

# Degree of the symmetric power X = Sym^n(C)
n = 15

print("Problem: Find the smallest and largest possible rank of the Neron-Severi group of X = Sym^15(C), where C is a genus 3 curve.")
print(f"Let g be the genus of C, g = {g}.")
print(f"Let n be the degree of the symmetric power, n = {n}.")

# Step 1: Use the relationship between the Neron-Severi ranks for a symmetric power
# and the Jacobian of the curve.
# For n > 2*g - 2, we have the formula: rho(Sym^n(C)) = rho(J(C)) + 1.
# Here, J(C) is the Jacobian of C.
condition_value = 2 * g - 2
print(f"\nWe check the condition n > 2*g - 2, which is {n} > {condition_value}. This is true.")
print("The formula rho(X) = rho(J(C)) + 1 applies.")

# Step 2: Determine the smallest possible rank for rho(J(C)).
# For a generic curve C, the rank of the Neron-Severi group of its Jacobian J(C) is minimal.
min_rho_j_c = 1
print(f"\nThe minimum possible rank for rho(J(C)) for a genus {g} curve is {min_rho_j_c}.")

# Step 3: Calculate the smallest rank for X.
smallest_rank = min_rho_j_c + 1
print(f"Therefore, the smallest possible rank of NS(X) is {min_rho_j_c} + 1 = {smallest_rank}.")

# Step 4: Determine the largest possible rank for rho(J(C)).
# The rank rho(J(C)) is bounded by g^2. This maximum is achieved for special curves
# with Complex Multiplication.
max_rho_j_c = g**2
print(f"\nThe maximum possible rank for rho(J(C)) for a genus {g} curve is g^2 = {g**2}, which is {max_rho_j_c}.")

# Step 5: Calculate the largest rank for X.
largest_rank = max_rho_j_c + 1
print(f"Therefore, the largest possible rank of NS(X) is {max_rho_j_c} + 1 = {largest_rank}.")
