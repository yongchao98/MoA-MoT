import fractions

def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    10-node barbell graph during uniform randomized gossiping.
    """
    # Total number of nodes in the graph
    total_nodes = 10

    # A barbell graph with 10 nodes consists of two K5 complete graphs
    # connected by a single edge.
    clique_size = 5

    # Let the bottleneck edge connect node 'u' and node 'v'.
    # Calculate the degree of the nodes connected by the bottleneck edge.
    # Each node is connected to (clique_size - 1) other nodes in its clique,
    # plus 1 connection via the bottleneck edge.
    degree_of_bottleneck_node = (clique_size - 1) + 1

    # Probability is calculated based on the two-step gossiping model:
    # 1. Select a node 'x' uniformly at random: P(x) = 1 / total_nodes
    # 2. 'x' selects a neighbor 'y' uniformly at random: P(y|x) = 1 / deg(x)

    # The bottleneck edge (u, v) is sampled if:
    #   - Node 'u' is picked AND 'u' gossips with 'v', OR
    #   - Node 'v' is picked AND 'v' gossips with 'u'.

    prob_u_picked = f"1 / {total_nodes}"
    prob_u_picks_v = f"1 / {degree_of_bottleneck_node}"

    prob_v_picked = f"1 / {total_nodes}"
    prob_v_picks_v = f"1 / {degree_of_bottleneck_node}"

    # Total Probability = P(u picked) * P(u picks v) + P(v picked) * P(v picks u)
    # The numbers are: 1/10, 1/5, 1/10, 1/5
    
    # Calculate the numerical result
    numerator = 2
    denominator = total_nodes * degree_of_bottleneck_node
    final_fraction = fractions.Fraction(numerator, denominator)


    print("Problem: Find the probability of sampling the bottleneck edge in a 10-node barbell graph via uniform gossiping.")
    print("\nStep 1: Define the graph structure.")
    print(f"A 10-node barbell graph has two complete graphs (K{clique_size}) connected by one 'bottleneck' edge.")

    print("\nStep 2: Determine the degrees of the nodes connected by the bottleneck edge.")
    print(f"A node on the bottleneck is connected to {clique_size - 1} nodes in its own clique and 1 node across the bridge.")
    print(f"So, its degree is ({clique_size - 1}) + 1 = {degree_of_bottleneck_node}.")

    print("\nStep 3: Calculate the probability based on the gossiping model.")
    print("The probability P(edge) = P(node 'u' is picked) * P('u' picks 'v') + P(node 'v' is picked) * P('v' picks 'u').")
    
    print("\nFinal Equation:")
    print(f"P(bottleneck) = (1 / {total_nodes}) * (1 / {degree_of_bottleneck_node}) + (1 / {total_nodes}) * (1 / {degree_of_bottleneck_node})")
    print(f"P(bottleneck) = (1/50) + (1/50)")
    print(f"P(bottleneck) = {numerator} / {denominator}")
    print(f"P(bottleneck) = {final_fraction.numerator} / {final_fraction.denominator}")
    
    final_decimal = float(final_fraction)
    print(f"P(bottleneck) = {final_decimal}")


solve_gossip_probability()
<<<0.04>>>