import sys
from fractions import Fraction

# Step 1 & 2: Define the graph parameters
total_nodes = 10
# A 10-node barbell graph has two cliques.
nodes_in_clique = total_nodes // 2

# Step 3: Calculate the degree of the nodes connected by the bottleneck.
# The degree of a node within a K_5 clique is 5-1 = 4.
degree_within_clique = nodes_in_clique - 1
# The bottleneck edge adds one more connection to these nodes.
degree_bottleneck_node = degree_within_clique + 1

# Step 4: Define the components of the probability formula.
# In uniform gossiping, any node is picked with probability 1/|V|.
prob_pick_node = Fraction(1, total_nodes)
# The selected node then picks a neighbor with probability 1/degree.
prob_contact_neighbor = Fraction(1, degree_bottleneck_node)

# Step 5: Calculate the final probability.
# The bottleneck edge (u,v) is sampled if u picks v OR v picks u.
# P = P(pick u)*P(u->v) + P(pick v)*P(v->u)
prob_u_picks_v = prob_pick_node * prob_contact_neighbor
prob_v_picks_u = prob_pick_node * prob_contact_neighbor
total_prob = prob_u_picks_v + prob_v_picks_u

# --- Output ---
print(f"For a barbell graph with {total_nodes} nodes, the structure is two {nodes_in_clique}-node complete graphs connected by one bottleneck edge.")
print(f"Let the bottleneck edge connect node 'u' and node 'v'.")
print(f"The total number of nodes |V| = {total_nodes}.")
print(f"The degree of the nodes 'u' and 'v' connected by the bottleneck is {degree_bottleneck_node}.")
print("\nThe probability of sampling the bottleneck edge in one step of uniform gossiping is calculated as:")
print("P = (Probabilty of picking u) * (Probability of u contacting v) + (Probabilty of picking v) * (Probability of v contacting u)")
print("P = (1/|V|) * (1/deg(u)) + (1/|V|) * (1/deg(v))")
print("\nSubstituting the values into the equation:")
print(f"P = (1/{total_nodes}) * (1/{degree_bottleneck_node}) + (1/{total_nodes}) * (1/{degree_bottleneck_node})")
print(f"P = {prob_u_picks_v} + {prob_v_picks_u}")
print(f"P = {total_prob}")
print(f"The final probability is {float(total_prob)}.")

# Suppress final answer line in notebook environments
if 'ipykernel' not in sys.modules:
    print(f"\n<<<{float(total_prob)}>>>")