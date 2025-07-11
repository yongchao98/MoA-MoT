import fractions

# Define the parameters of the barbell graph
total_nodes = 10
# The graph consists of two equal-sized cliques
clique_size = total_nodes // 2

# Calculate the degree of the two nodes connected by the bottleneck edge.
# Each of these nodes is connected to (clique_size - 1) nodes within its clique,
# plus 1 connection via the bottleneck edge.
degree_of_bottleneck_node = (clique_size - 1) + 1

# In randomized uniform gossiping, the probability of sampling an edge (A, B) is:
# P(A is chosen) * P(A chooses B) + P(B is chosen) * P(B chooses A)
# which is (1/N) * (1/deg(A)) + (1/N) * (1/deg(B))

# Calculate the probability
prob_pick_node = fractions.Fraction(1, total_nodes)
prob_node_chooses_neighbor = fractions.Fraction(1, degree_of_bottleneck_node)

prob_one_way = prob_pick_node * prob_node_chooses_neighbor
total_prob = 2 * prob_one_way

# Print the step-by-step derivation of the formula
print("The probability 'P' of sampling the bottleneck edge is the sum of two cases:")
print("1. Picking the first bottleneck node and it choosing the second.")
print("2. Picking the second bottleneck node and it choosing the first.")
print("\nP = P(pick node A) * P(A chooses B) + P(pick node B) * P(B chooses A)")
print(f"P = (1/{total_nodes}) * (1/{degree_of_bottleneck_node}) + (1/{total_nodes}) * (1/{degree_of_bottleneck_node})")
print(f"P = {prob_one_way} + {prob_one_way}")
print(f"P = {total_prob}")
print(f"\nThe final probability is {float(total_prob)}")