import sys

def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a symmetric
    10-node barbell graph during randomized uniform gossiping.
    """
    # Total number of nodes in the graph
    total_nodes = 10

    # For a symmetric barbell graph, the nodes are split into two equal-sized cliques.
    nodes_per_clique = int(total_nodes / 2)

    # The two cliques are connected by a single bottleneck edge.
    # Let's find the degree of a node connected to this edge.
    # It is connected to all other nodes in its clique (nodes_per_clique - 1)
    # plus the one node across the bottleneck.
    degree_of_bottleneck_node = (nodes_per_clique - 1) + 1

    # The probability of selecting any specific node in the graph is 1 / total_nodes.
    prob_select_node = 1 / total_nodes

    # The probability of a selected bottleneck node choosing its partner across the bridge
    # is 1 / its degree.
    prob_select_neighbor = 1 / degree_of_bottleneck_node

    # The sampling of the bottleneck edge can happen in two ways:
    # 1. Node u is picked, and it picks neighbor v.
    prob_case1 = prob_select_node * prob_select_neighbor
    # 2. Node v is picked, and it picks neighbor u.
    prob_case2 = prob_select_node * prob_select_neighbor

    # The total probability is the sum of these two mutually exclusive cases.
    total_probability = prob_case1 + prob_case2

    # The final equation is P(total) = P(u->v) + P(v->u)
    # P(u->v) = P(u is picked) * P(v is picked | u was picked)
    # P(v->u) = P(v is picked) * P(u is picked | v was picked)
    # The numbers are: 1/total_nodes, 1/degree_of_bottleneck_node for each case.

    # Redirecting print to sys.stdout as per platform convention.
    original_stdout = sys.stdout
    try:
        # Final outputting of the equation with each number.
        print("The final probability is calculated from the sum of two scenarios:")
        print("1. Picking one bottleneck node, then its partner: (1/10) * (1/5)")
        print("2. Picking the other bottleneck node, then its partner: (1/10) * (1/5)")
        print("\nThe final equation is:")
        print(f"({prob_select_node:.2f} * {prob_select_neighbor:.2f}) + ({prob_select_node:.2f} * {prob_select_neighbor:.2f}) = {total_probability:.2f}")

        # Printing the fractional representation for clarity
        print(f"1/50 + 1/50 = 2/50 = 1/25")

    finally:
        sys.stdout = original_stdout


solve_gossip_probability()