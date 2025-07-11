# Total number of nodes in the barbell graph
total_nodes = 10

# The barbell graph consists of two cliques.
# Size of each clique (m) is total_nodes / 2.
clique_size = total_nodes / 2

# In a barbell graph, the two cliques are connected by a single "bottleneck" edge.
# Let the nodes at the ends of this bottleneck edge be u_A and u_B.

# The degree of a node at the end of the bottleneck edge is calculated as:
# (connections inside its clique) + (the bottleneck connection)
# Inside a clique of size m, a node is connected to m-1 other nodes.
# So, degree = (m - 1) + 1 = m.
degree_bottleneck_node = int(clique_size)

# The gossiping process involves:
# 1. Picking a node `u` uniformly at random from all N nodes. P(pick u) = 1/N.
# 2. The chosen node `u` picks a neighbor `v` uniformly at random. P(u picks v) = 1/degree(u).

# We want to find the probability of sampling the bottleneck edge (u_A, u_B).
# This can happen in two ways:
# 1. Node u_A is picked, and it picks its neighbor u_B.
prob_a_picks_b = (1 / total_nodes) * (1 / degree_bottleneck_node)

# 2. Node u_B is picked, and it picks its neighbor u_A.
prob_b_picks_a = (1 / total_nodes) * (1 / degree_bottleneck_node)

# The total probability is the sum of these two mutually exclusive events.
total_probability = prob_a_picks_b + prob_b_picks_a

# Print the step-by-step calculation
print("A barbell graph with 10 nodes has two cliques of 5 nodes each.")
print(f"The two nodes at the ends of the connecting bottleneck edge each have a degree of {degree_bottleneck_node}.")
print("\nThe probability of sampling the bottleneck edge is the sum of two scenarios:")
print("1. Picking the first bottleneck node AND that node picking the second.")
print("2. Picking the second bottleneck node AND that node picking the first.")
print("\nFinal Calculation:")
print(f"P(bottleneck) = P(pick u_A) * P(u_A picks u_B) + P(pick u_B) * P(u_B picks u_A)")
print(f"P(bottleneck) = (1 / {total_nodes}) * (1 / {degree_bottleneck_node}) + (1 / {total_nodes}) * (1 / {degree_bottleneck_node})")
print(f"P(bottleneck) = {prob_a_picks_b} + {prob_b_picks_a}")
print(f"P(bottleneck) = {total_probability}")
print(f"As a fraction, this is 2/50 or 1/25.")
