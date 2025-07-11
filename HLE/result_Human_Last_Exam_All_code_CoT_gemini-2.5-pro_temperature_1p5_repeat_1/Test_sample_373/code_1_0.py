def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes during uniform gossiping.
    """
    # Total number of nodes in the barbell graph
    total_nodes = 10

    # A 10-node barbell graph consists of two 5-node complete graphs (cliques).
    nodes_in_clique = total_nodes // 2

    # The bottleneck edge connects one node from each clique. Let's call them u and v.
    # We need to calculate their degrees.
    # The degree of a node in a K_n graph is n-1.
    # The bottleneck edge adds 1 to the degree of u and v.
    degree_u = (nodes_in_clique - 1) + 1
    degree_v = (nodes_in_clique - 1) + 1

    # The probability of sampling an edge (u, v) in uniform gossiping is:
    # P(sample) = P(choose u) * P(u chooses v) + P(choose v) * P(v chooses u)
    # P(choose node) = 1 / total_nodes
    # P(node chooses neighbor) = 1 / degree(node)

    prob_u_chooses_v = (1 / total_nodes) * (1 / degree_u)
    prob_v_chooses_u = (1 / total_nodes) * (1 / degree_v)

    total_probability = prob_u_chooses_v + prob_v_chooses_u

    # Output the final equation with all the numbers, as requested.
    print(f"The barbell graph has N = {total_nodes} nodes.")
    print(f"It consists of two K_{nodes_in_clique} complete graphs.")
    print(f"The two nodes on the bottleneck edge, u and v, have degrees:")
    print(f"deg(u) = {degree_u}, deg(v) = {degree_v}")
    print("\nThe probability of sampling the bottleneck edge (u, v) is:")
    print(f"P(sample) = P(choose u) * P(u chooses v) + P(choose v) * P(v chooses u)")
    print(f"P(sample) = (1/{total_nodes}) * (1/{degree_u}) + (1/{total_nodes}) * (1/{degree_v})")
    print(f"P(sample) = ({prob_u_chooses_v}) + ({prob_v_chooses_u})")
    print(f"P(sample) = {total_probability}")
    print("\nAs a fraction, this is 1/50 + 1/50 = 2/50 = 1/25.")

solve_gossip_probability()