import fractions

def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes during randomized uniform gossiping.
    """
    # 1. Define the graph parameters
    # Total number of nodes in the barbell graph
    total_nodes_n = 10
    # Size of each complete graph (clique)
    clique_size_m = 5

    # 2. Calculate the degree of the bridge nodes (u and v)
    # The degree is the number of nodes in its clique (m-1) plus the bridge connection (1)
    degree_bridge_node = (clique_size_m - 1) + 1

    # 3. Calculate the probability of sampling the edge in each direction
    # P(u -> v) = P(picking u) * P(picking v | u) = (1/n) * (1/deg(u))
    prob_u_to_v = fractions.Fraction(1, total_nodes_n) * fractions.Fraction(1, degree_bridge_node)
    
    # P(v -> u) = P(picking v) * P(picking u | v) = (1/n) * (1/deg(v))
    prob_v_to_u = fractions.Fraction(1, total_nodes_n) * fractions.Fraction(1, degree_bridge_node)

    # 4. The total probability is the sum of the probabilities of traversing the edge in either direction
    total_probability = prob_u_to_v + prob_v_to_u

    # 5. Print the explanation and the final equation
    print("A 10-node barbell graph consists of two 5-node cliques connected by a single 'bottleneck' edge.")
    print(f"Let the bridge nodes be 'u' and 'v'.")
    print(f"Total nodes n = {total_nodes_n}")
    print(f"Degree of bridge node u, deg(u) = (clique_size - 1) + 1 = ({clique_size_m} - 1) + 1 = {degree_bridge_node}")
    print(f"Degree of bridge node v, deg(v) = (clique_size - 1) + 1 = ({clique_size_m} - 1) + 1 = {degree_bridge_node}")
    print("\nThe probability of sampling the bottleneck edge P(edge_uv) is the sum of probabilities of gossiping in either direction:")
    print("P(edge_uv) = P(u -> v) + P(v -> u)")
    print("P(edge_uv) = (P(pick u) * P(pick v | u)) + (P(pick v) * P(pick u | v))")
    print(f"P(edge_uv) = (1/n * 1/deg(u)) + (1/n * 1/deg(v))")
    print(f"P(edge_uv) = (1/{total_nodes_n} * 1/{degree_bridge_node}) + (1/{total_nodes_n} * 1/{degree_bridge_node})")
    print(f"P(edge_uv) = {prob_u_to_v} + {prob_v_to_u}")
    print(f"P(edge_uv) = {total_probability}")
    print(f"P(edge_uv) = {float(total_probability)}")

solve_gossip_probability()