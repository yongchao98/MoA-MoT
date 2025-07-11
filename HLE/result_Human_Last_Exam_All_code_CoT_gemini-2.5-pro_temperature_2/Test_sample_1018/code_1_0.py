import math

# Genera of the two surfaces
g1 = 31
g2 = 17

# Dimensions of the surfaces
n1 = 2
n2 = 2

# Step 1: Compute the simplicial volume of the first surface, Sigma_31
# Formula: ||Sigma_g|| = 4g - 4 for g >= 2
sv1 = 4 * g1 - 4

# Step 2: Compute the simplicial volume of the second surface, Sigma_17
# Formula: ||Sigma_g|| = 4g - 4 for g >= 2
sv2 = 4 * g2 - 4

# Step 3: Compute the binomial coefficient for the product formula
# C(n1, n2) = (n1+n2)! / (n1! * n2!)
binomial_coeff = math.comb(n1 + n2, n1)

# Step 4: Compute the final simplicial volume of the product space
# ||M1 x M2|| = C(n1, n2) * ||M1|| * ||M2||
final_sv = binomial_coeff * sv1 * sv2

# Print the final equation with all the computed components
print(f"The simplicial volume of Sigma_{g1} x Sigma_{g2} is calculated as:")
print(f"||Sigma_{g1} x Sigma_{g2}|| = C({n1}, {n2}) * ||Sigma_{g1}|| * ||Sigma_{g2}||")
print(f"Which evaluates to:")
print(f"result = {binomial_coeff} * {sv1} * {sv2} = {final_sv}")