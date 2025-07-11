import sympy

# Define the symbols
K, K1, K2, K3, K4 = sympy.symbols('K K1 K2 K3 K4')

# Values derived from the analysis
k_val = -1
k1_val = 4
k2_val = 2
k3_val = -1
k4_val = 2

# Calculate the product
product = k_val * k1_val * k2_val * k3_val * k4_val

# Print the equation and the result
print(f"The product we want to calculate is K * K1 * K2 * K3 * K4.")
print(f"Based on the analysis, we found the values:")
print(f"K = {k_val}")
print(f"K1 = {k1_val}")
print(f"K2 = {k2_val}")
print(f"K3 = {k3_val}")
print(f"K4 = {k4_val}")
print("\nSubstituting these values into the product:")
print(f"Product = ({k_val}) * ({k1_val}) * ({k2_val}) * ({k3_val}) * ({k4_val}) = {product}")