# Define the parameters of the barbell graph
total_nodes_N = 10
# For a 10-node barbell, the most symmetric structure is two 5-node cliques
clique_size_k = 5

# The barbell graph is formed by connecting two complete graphs (cliques) K_k
# with a single "bottleneck" edge. Let the nodes of this edge be u and v.

# Calculate the degree of the nodes connected by the bottleneck edge.
# The degree is the number of connections within the clique (k-1) plus the
# one connection through the bottleneck edge.
degree_u = (clique_size_k - 1) + 1
degree_v = degree_u # By symmetry

# In randomized uniform gossiping, a node is chosen uniformly at random (with probability 1/N),
# and then one of its neighbors is chosen uniformly at random (with probability 1/degree).

# The probability of sampling the bottleneck edge (u, v) is the sum of two cases:
# 1. Node u is chosen, and it picks neighbor v.
# 2. Node v is chosen, and it picks neighbor u.
# The formula is: P = (1/N) * (1/d(u)) + (1/N) * (1/d(v))

# Calculate the final probability
prob_u_picks_v = (1 / total_nodes_N) * (1 / degree_u)
prob_v_picks_u = (1 / total_nodes_N) * (1 / degree_v)
total_prob = prob_u_picks_v + prob_v_picks_u

print("To find the probability of sampling the bottleneck edge in a 10-node barbell graph:")
print(f"1. Define the graph: Two {clique_size_k}-node cliques connected by one edge.")
print(f"   - Total nodes N = {total_nodes_N}")
print(f"   - Let the bottleneck edge be (u, v).")
print(f"2. Calculate the degree of nodes u and v:")
print(f"   - Degree d(u) = d(v) = (connections in clique) + (bottleneck) = ({clique_size_k}-1) + 1 = {degree_u}")
print("3. Calculate the probability using the gossiping formula: P = (1/N)*(1/d(u)) + (1/N)*(1/d(v))")
print("\nFinal Equation:")
print(f"P = (1/{total_nodes_N}) * (1/{degree_u}) + (1/{total_nodes_N}) * (1/{degree_v})")
print(f"P = {prob_u_picks_v} + {prob_v_picks_u}")
print(f"P = {total_prob}")
