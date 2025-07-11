import math

def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a barbell graph
    during a randomized uniform gossiping step.
    """
    # Problem parameters
    total_nodes = 10

    # Explanation of the plan
    print("This script calculates the probability of sampling the bottleneck edge in a 10-node barbell graph.")
    print("-" * 50)

    # Step 1: Define the graph's structure.
    print("\nStep 1: Define the Graph Structure")
    clique_size = total_nodes // 2
    print(f"A barbell graph with {total_nodes} nodes consists of two complete graphs (K_{clique_size}) connected by a single bottleneck edge.")

    # Step 2: Determine the degrees of the nodes connected to the bottleneck.
    print("\nStep 2: Determine Node Degrees")
    # In a complete graph K_k, each node has a degree of k-1.
    degree_in_clique = clique_size - 1
    # The two nodes at the ends of the bottleneck edge have one additional connection.
    degree_bridge_node = degree_in_clique + 1
    print(f"The two nodes connected by the bottleneck each belong to a K_{clique_size}.")
    print(f"Their degree within their clique is {clique_size} - 1 = {degree_in_clique}.")
    print(f"Adding the bottleneck connection, their total degree is {degree_in_clique} + 1 = {degree_bridge_node}.")

    # Step 3: Calculate the probability.
    print("\nStep 3: Calculate the Probability")
    print("The bottleneck edge can be selected in two ways during a gossip step:")
    print("  1. Start at the bridge node in the first clique and select its neighbor in the second clique.")
    print("  2. Start at the bridge node in the second clique and select its neighbor in the first clique.")

    prob_select_node = 1 / total_nodes
    prob_select_neighbor = 1 / degree_bridge_node
    
    # Total probability is the sum of the two mutually exclusive scenarios.
    total_prob = (prob_select_node * prob_select_neighbor) + (prob_select_node * prob_select_neighbor)
    
    total_prob_num = 2
    total_prob_den = total_nodes * degree_bridge_node

    common_divisor = math.gcd(total_prob_num, total_prob_den)
    simple_num = total_prob_num // common_divisor
    simple_den = total_prob_den // common_divisor

    print("\nThe probability equation, showing each number, is as follows:")
    print(f"P(bottleneck) = [P(select node 1) * P(select neighbor 2)] + [P(select node 2) * P(select neighbor 1)]")
    print(f"              = [ (1 / {total_nodes}) * (1 / {degree_bridge_node}) ] + [ (1 / {total_nodes}) * (1 / {degree_bridge_node}) ]")
    print(f"              = (1 / {total_nodes * degree_bridge_node}) + (1 / {total_nodes * degree_bridge_node})")
    print(f"              = {total_prob_num} / {total_prob_den}")
    print(f"              = {simple_num} / {simple_den}")
    print(f"              = {total_prob}")

solve_gossip_probability()
<<<0.04>>>