import math
from fractions import Fraction

# Plan:
# 1. Define graph parameters for a 10-node barbell graph.
# 2. Calculate the degree of the two nodes connected by the bottleneck edge.
# 3. Use the formula for uniform gossiping probability to find the probability of
#    traversing the bottleneck edge in either direction.
# 4. Sum the probabilities and present the final equation and result.

# 1. Graph Parameters
# Total number of nodes in the graph
total_nodes = 10
# A symmetric barbell graph with 10 nodes consists of two complete graphs
# of size m, connected by a single bottleneck edge.
nodes_per_clique = total_nodes // 2

# 2. Node Degree Calculation
# The two nodes connected by the bottleneck edge are 'u' and 'v'.
# The degree of a bridge node 'u' is its connections within its own clique (m-1)
# plus the single connection across the bridge.
degree_bridge_node = (nodes_per_clique - 1) + 1

# 3. Probability Calculation
# The probability of traversing the directed edge (u, v) is P(u) * P(v|u).
# P(u -> v) = (1/N) * (1/degree(u))
prob_u_to_v = Fraction(1, total_nodes * degree_bridge_node)

# By symmetry, the probability of traversing (v, u) is the same.
prob_v_to_u = Fraction(1, total_nodes * degree_bridge_node)

# The total probability is the sum, as the bottleneck edge can be
# sampled by traversing it in either direction.
total_prob = prob_u_to_v + prob_v_to_u

# 4. Output Results
print(f"The probability of sampling the bottleneck edge is calculated as follows:")
print(f"The graph has N = {total_nodes} nodes, composed of two K_{nodes_per_clique} cliques.")
print(f"Let the bottleneck edge connect nodes 'u' and 'v'.")
print(f"The degree of a bridge node is: Degree(u) = Degree(v) = ({nodes_per_clique} - 1) + 1 = {degree_bridge_node}")
print("")
print("The probability of sampling the edge is the sum of probabilities of traversing it in either direction:")
print("P(bottleneck) = P(u -> v) + P(v -> u)")
print("P(bottleneck) = [P(pick u) * P(pick v | u)] + [P(pick v) * P(pick u | v)]")
print("P(bottleneck) = [ (1/N) * (1/Degree(u)) ] + [ (1/N) * (1/Degree(v)) ]")
print("")
print("Substituting the values into the equation:")
print(f"P(bottleneck) = (1 / {total_nodes} * 1 / {degree_bridge_node}) + (1 / {total_nodes} * 1 / {degree_bridge_node})")
denominator = total_nodes * degree_bridge_node
print(f"P(bottleneck) = (1 / {denominator}) + (1 / {denominator})")
print(f"P(bottleneck) = {total_prob.numerator} / {total_prob.denominator}")
print(f"P(bottleneck) = {float(total_prob)}")