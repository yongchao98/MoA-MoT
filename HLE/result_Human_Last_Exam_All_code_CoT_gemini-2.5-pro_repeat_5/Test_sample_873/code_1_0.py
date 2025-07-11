# Step 1: Define the parameters for the graph G = Km - E(C5).
# The Shannon capacity of G is determined by the removed 5-cycle (C5).
# |V(C5)| is the number of vertices in the 5-cycle.
num_vertices_C5 = 5
# α(C5) is the independence number of the 5-cycle (largest set of non-adjacent vertices).
alpha_C5 = 2

# Step 2: Calculate the Shannon capacity of G using Haemers' theorem.
# c(G) = |V(C5)| / α(C5)
shannon_capacity_G = num_vertices_C5 / alpha_C5

# Step 3: Define the parameters for the graph H = Kn - E(C4).
# The Shannon capacity of H is determined by the removed 4-cycle (C4).
# |V(C4)| is the number of vertices in the 4-cycle.
num_vertices_C4 = 4
# α(C4) is the independence number of the 4-cycle.
alpha_C4 = 2

# Step 4: Calculate the Shannon capacity of H.
# c(H) = |V(C4)| / α(C4)
shannon_capacity_H = num_vertices_C4 / alpha_C4

# Step 5: The Shannon capacity of the strong product G⊠H is the product of their individual capacities.
# c(G⊠H) = c(G) * c(H)
final_capacity = shannon_capacity_G * shannon_capacity_H

# Step 6: Print the final equation and the result.
# The problem asks to output each number in the final equation.
print(f"The Shannon capacity of G is |V(C5)| / α(C5) = {num_vertices_C5} / {alpha_C5} = {shannon_capacity_G}")
print(f"The Shannon capacity of H is |V(C4)| / α(C4) = {num_vertices_C4} / {alpha_C4} = {shannon_capacity_H}")
print(f"The Shannon capacity of G⊠H is c(G) * c(H).")
print(f"The calculation is: ({num_vertices_C5} / {alpha_C5}) * ({num_vertices_C4} / {alpha_C4}) = {final_capacity}")