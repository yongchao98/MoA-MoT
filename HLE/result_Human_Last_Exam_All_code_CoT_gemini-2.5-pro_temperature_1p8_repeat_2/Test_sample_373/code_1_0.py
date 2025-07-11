import fractions

# Step 1: Define the parameters of the barbell graph.
total_nodes = 10
# A standard barbell graph splits the nodes into two equal-sized cliques.
clique_size = total_nodes // 2

# Step 2: Calculate the degree of the two nodes connected by the bottleneck edge.
# A node in a complete graph (clique) of size n has degree n-1 within that clique.
degree_internal = clique_size - 1
# The nodes connected to the bottleneck have their internal connections plus the bridge.
degree_bottleneck_node = degree_internal + 1

# Step 3: In uniform gossiping, the probability of selecting a neighbor is 1/degree.
# Let P(A) be the probability of the first node choosing the second.
prob_A = fractions.Fraction(1, degree_bottleneck_node)
# Let P(B) be the probability of the second node choosing the first.
prob_B = fractions.Fraction(1, degree_bottleneck_node)

# Step 4: The bottleneck edge is sampled if A or B happens.
# Since the choices are independent, P(A or B) = P(A) + P(B) - P(A) * P(B).
total_prob = prob_A + prob_B - (prob_A * prob_B)

# Print the explanation and the final result.
print(f"A barbell graph with {total_nodes} nodes has two cliques of size {clique_size}.")
print(f"The nodes connecting the cliques each have a degree of ({clique_size} - 1) + 1 = {degree_bottleneck_node}.")
print("\nThe probability of sampling the bottleneck edge is the probability that either node connected by it chooses the other.")
print("This is calculated as P(A) + P(B) - P(A)*P(B), where A and B are the events of each node choosing the other.")
print("\nFinal Equation:")
# Output each number in the final equation as requested.
print(f"P(bottleneck sampled) = {prob_A} + {prob_B} - ({prob_A} * {prob_B}) = {total_prob}")
print(f"\nThe decimal value of this probability is {float(total_prob)}.")
