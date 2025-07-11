def calculate_bottleneck_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a barbell graph
    with 10 nodes using randomized uniform gossiping.
    """
    # 1. Define graph parameters
    total_nodes = 10
    # A barbell graph has two equal-sized cliques
    clique_size = total_nodes // 2

    # 2. Calculate the degree of the nodes connected by the bottleneck
    # A node at the bottleneck is connected to (clique_size - 1) nodes in its
    # own clique, plus 1 node in the other clique via the bottleneck edge.
    degree_u = (clique_size - 1) + 1
    degree_v = (clique_size - 1) + 1

    # 3. Calculate the probability
    # The probability of sampling edge (u, v) is:
    # P = P(pick u) * P(u picks v) + P(pick v) * P(v picks u)
    # P = (1/total_nodes) * (1/degree_u) + (1/total_nodes) * (1/degree_v)
    prob_pick_u_then_v = (1 / total_nodes) * (1 / degree_u)
    prob_pick_v_then_u = (1 / total_nodes) * (1 / degree_v)
    total_prob = prob_pick_u_then_v + prob_pick_v_then_u

    # 4. Print the final equation with the numbers and the result
    print("The probability of sampling the bottleneck edge (u,v) is calculated as:")
    print("P = (1/Total_Nodes) * (1/Degree(u)) + (1/Total_Nodes) * (1/Degree(v))")
    print("\nPlugging in the values for a 10-node barbell graph:")
    print(f"P = (1/{total_nodes}) * (1/{degree_u}) + (1/{total_nodes}) * (1/{degree_v})")
    print(f"\nThis simplifies to:")
    # We use integers for the denominator to show the fractional form clearly
    final_denominator_u = total_nodes * degree_u
    final_denominator_v = total_nodes * degree_v
    print(f"P = 1/{final_denominator_u} + 1/{final_denominator_v}")
    print("\nThe final probability is:")
    print(total_prob)


if __name__ == "__main__":
    calculate_bottleneck_probability()
<<<0.04>>>