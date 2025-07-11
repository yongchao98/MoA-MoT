import math

# Step 1: Define the degree 'k' of the central node.
# The graph has 64 nodes, and the central node is connected to all other nodes.
k = 64 - 1

# Step 2: Define the number of edges 'E' between the neighbors.

# Edges among the 3 "blue" neighbors (they form a K3 clique).
edges_blue = 3

# Edges within the 3 "green" K4-cliques. Each K4 has 6 edges.
edges_green = 3 * 6

# The total number of edges among neighbors in the central 16-node structure.
edges_in_central_structure = edges_blue + edges_green

# Edges in the 3 "red" structures. Each is a copy of the central one.
edges_red = 3 * edges_in_central_structure

# Total edges 'E' between all neighbors of the central node.
E = edges_in_central_structure + edges_red

# Step 3: Calculate the clustering coefficient 'C'.
if k < 2:
    C = 0.0
else:
    C_numerator = 2 * E
    C_denominator = k * (k - 1)
    C = C_numerator / C_denominator

# Find the simplified fraction
common_divisor = math.gcd(C_numerator, C_denominator)
simple_num = C_numerator // common_divisor
simple_den = C_denominator // common_divisor


# Step 4: Print the explanation, the formula with numbers, and the final result.
print("Calculation of the Clustering Coefficient")
print("---------------------------------------")
print("Formula: C = (2 * E) / (k * (k - 1))")
print(f"\n1. The degree 'k' of the central node is the total number of nodes minus one:")
print(f"   k = 64 - 1 = {k}")

print(f"\n2. The number of edges 'E' between these {k} neighbors is calculated as follows:")
print(f"   - Edges from the 3 blue neighbors (a triangle): {edges_blue}")
print(f"   - Edges from the 3 green groups (3 * 6): {edges_green}")
print(f"   - Edges from the 3 red structures (3 * ({edges_blue} + {edges_green})): {edges_red}")
print(f"   Total E = {edges_blue} + {edges_green} + {edges_red} = {E}")

print(f"\n3. Plugging these values into the formula:")
print(f"   C = (2 * {E}) / ({k} * ({k} - 1))")
print(f"   C = {C_numerator} / ({k} * {k-1})")
print(f"   C = {C_numerator} / {C_denominator}")

print(f"\nThe simplified fraction is {simple_num}/{simple_den}.")
print(f"The final clustering coefficient is approximately: {C}")
