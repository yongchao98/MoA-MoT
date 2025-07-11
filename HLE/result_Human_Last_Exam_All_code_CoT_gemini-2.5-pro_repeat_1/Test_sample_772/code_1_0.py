import math

# Step 1: Define the parameters for the group G = SO_3(R)
# The dimension of the group G
d = 3.0

# Step 2: Define the dimensions of possible proper Lie subgroups H of G
# Subgroups can be finite (like {e}) with dimension 0, or isomorphic to SO(2) with dimension 1.
d_H_options = [0.0, 1.0]

# Step 3: Define a function to calculate the growth exponent E
def calculate_growth_exponent(d_H, d):
    """
    Calculates the growth exponent E for a neighborhood of a subgroup H.
    The model assumes ballistic growth along H and diffusive growth orthogonally.
    The measure of the n-th product set X^n grows as n^E.
    E = d_H + (d - d_H) / 2
    """
    return d_H + (d - d_H) / 2.0

# Step 4: Calculate the exponent for each possible subgroup dimension
print("Analyzing the growth of product sets in G = SO_3(R).")
print(f"The dimension of the group is d = {int(d)}.")
print(f"The possible dimensions for proper subgroups H are d_H = {[int(h) for h in d_H_options]}.")
print("-" * 20)
print("The growth of measure for a neighborhood X of H follows mu(X^n) ~ n^E * mu(X).")
print("The growth exponent E is calculated as E = d_H + (d - d_H) / 2.")
print("-" * 20)

exponents = {}
for d_H in d_H_options:
    exponent = calculate_growth_exponent(d_H, d)
    exponents[d_H] = exponent
    print(f"For a subgroup H with dimension d_H = {int(d_H)}:")
    print(f"  The exponent E = {d_H} + ({d} - {d_H}) / 2 = {exponent}")

# Step 5: Find the minimum exponent, which corresponds to the slowest growth (worst case)
min_exponent = min(exponents.values())

# Step 6: Calculate alpha based on the minimum exponent
# The relation n(N) is determined by the slowest growth: n(N)^E_min * (1/N) ~ 1
# This gives n(N) ~ N^(1/E_min), so alpha = 1 / E_min.
alpha = 1.0 / min_exponent

print("-" * 20)
print(f"The slowest growth occurs for the minimum exponent, which is min(E) = {min_exponent}.")
print("This corresponds to a set that is a small ball-like neighborhood (d_H = 0).")
print("The scaling relation is n(N) ~ N^alpha, where alpha is given by the equation: alpha = 1 / min(E).")
print(f"The final calculation is: alpha = 1 / {min_exponent}")
print(f"The resulting value for alpha is: {alpha}")
