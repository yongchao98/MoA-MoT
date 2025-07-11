import fractions

# Graph parameters
total_nodes = 10
clique_size = 5

# The bottleneck edge connects two nodes, one from each clique.
# Let's calculate the degree of these two specific nodes.
# Degree = (connections within their clique) + (the bottleneck connection)
degree_of_bottleneck_node = (clique_size - 1) + 1

# --- Probability Calculation ---
# The gossiping process:
# 1. A node is selected uniformly at random from all nodes.
prob_select_any_node = fractions.Fraction(1, total_nodes)

# 2. The selected node chooses a neighbor uniformly at random.
# We only care if one of the two bottleneck nodes is chosen.
prob_bottleneck_node_selects_bridge = fractions.Fraction(1, degree_of_bottleneck_node)

# The bottleneck edge can be sampled in two symmetric ways:
# - Pick the first bottleneck node AND that node picks the bridge.
# - Pick the second bottleneck node AND that node picks thebridge.
prob_path1 = prob_select_any_node * prob_bottleneck_node_selects_bridge
prob_path2 = prob_select_any_node * prob_bottleneck_node_selects_bridge

total_probability = prob_path1 + prob_path2

# --- Output the results ---
print("Calculating the probability of sampling the bottleneck edge in a 10-node barbell graph.")
print(f"A 10-node barbell graph has two cliques of {clique_size} nodes.")
print(f"The degree of each of the 2 nodes connected by the bottleneck edge is {degree_of_bottleneck_node}.")
print("\nThe probability is the sum of two scenarios:")
print("1. Picking the first bottleneck node (prob 1/10) AND it choosing its bridge neighbor (prob 1/5).")
print("2. Picking the second bottleneck node (prob 1/10) AND it choosing its bridge neighbor (prob 1/5).")

print("\nFinal Equation:")
# The problem asks to output each number in the final equation.
# P = (prob_select_node_A * prob_A_selects_B) + (prob_select_node_B * prob_B_selects_A)
# prob_select_node_A and prob_select_node_B are both prob_select_any_node.
# prob_A_selects_B and prob_B_selects_A are both prob_bottleneck_node_selects_bridge.
print(f"P(bottleneck) = ({prob_select_any_node.numerator}/{prob_select_any_node.denominator} * "
      f"{prob_bottleneck_node_selects_bridge.numerator}/{prob_bottleneck_node_selects_bridge.denominator}) + "
      f"({prob_select_any_node.numerator}/{prob_select_any_node.denominator} * "
      f"{prob_bottleneck_node_selects_bridge.numerator}/{prob_bottleneck_node_selects_bridge.denominator})")

print(f"P(bottleneck) = {prob_path1.numerator}/{prob_path1.denominator} + "
      f"{prob_path2.numerator}/{prob_path2.denominator}")

print(f"P(bottleneck) = {total_probability.numerator}/{total_probability.denominator}")
print(f"\nThe probability as a decimal is: {float(total_probability)}")

<<<0.04>>>