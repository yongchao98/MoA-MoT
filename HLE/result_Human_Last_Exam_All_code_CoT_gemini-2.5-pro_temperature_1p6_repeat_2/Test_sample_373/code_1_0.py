import fractions

def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a 10-node
    barbell graph using randomized uniform gossiping.
    """
    # Step 1: Define the graph parameters.
    # A barbell graph with 10 nodes is composed of two 5-node cliques (K5)
    # connected by a single bottleneck edge.
    total_nodes = 10
    clique_size = 5

    # Step 2: Determine the degrees of the two nodes connected by the bottleneck edge.
    # Let the bottleneck edge connect node 'u' and node 'v'.
    # The degree of a node is its number of connections.
    # Each node is connected to (clique_size - 1) nodes within its own clique,
    # plus 1 connection via the bottleneck edge.
    degree_of_bridge_node = (clique_size - 1) + 1

    # Step 3: Formulate and calculate the probability.
    # The uniform gossiping process is:
    # 1. Pick a node 'i' uniformly at random from N total nodes (Probability = 1/N).
    # 2. Pick one of its neighbors 'j' uniformly at random (Probability = 1/degree(i)).
    #
    # The bottleneck edge (u, v) is sampled if:
    # - Node 'u' is picked AND 'u' chooses 'v'.
    # OR
    # - Node 'v' is picked AND 'v' chooses 'u'.
    #
    # The probability is P = P(pick u)*P(u chooses v) + P(pick v)*P(v chooses u)
    
    # P(pick a specific node)
    prob_pick_node = fractions.Fraction(1, total_nodes)
    
    # P(a bridge node chooses its counterpart across the bridge)
    prob_choose_across_bridge = fractions.Fraction(1, degree_of_bridge_node)
    
    # Probability of one side of the event (e.g., u picks v)
    prob_one_direction = prob_pick_node * prob_choose_across_bridge
    
    # Total probability is the sum of both directions
    total_prob = prob_one_direction + prob_one_direction

    # Step 4: Print the detailed explanation and the final answer.
    print("This script calculates the probability of sampling the bottleneck edge in a 10-node barbell graph.")
    print("\n--- Problem Setup ---")
    print(f"Total number of nodes (N): {total_nodes}")
    print("The graph consists of two 5-node cliques connected by one bottleneck edge.")
    print("Let the bottleneck edge connect node 'u' and node 'v'.")

    print("\n--- Degree Calculation ---")
    print("The degree of a node is its number of connections.")
    print(f"Degree of u (d_u) = (connections in clique) + (bottleneck) = {clique_size - 1} + 1 = {degree_of_bridge_node}")
    print(f"Degree of v (d_v) = (connections in clique) + (bottleneck) = {clique_size - 1} + 1 = {degree_of_bridge_node}")

    print("\n--- Probability Calculation ---")
    print("The probability 'P' is the sum of two scenarios:")
    print("1. Node 'u' is picked, and it chooses 'v'.")
    print("2. Node 'v' is picked, and it chooses 'u'.")
    print("\nThe formula is:")
    print("P = [P(pick u) * P(u chooses v)] + [P(pick v) * P(v chooses u)]")
    print("P = [ (1/N) * (1/d_u) ] + [ (1/N) * (1/d_v) ]")
    
    print("\nPlugging in the numbers:")
    print(f"P = [ (1/{total_nodes}) * (1/{degree_of_bridge_node}) ] + [ (1/{total_nodes}) * (1/{degree_of_bridge_node}) ]")
    print(f"P = [ {prob_one_direction} ] + [ {prob_one_direction} ]")
    print(f"P = {total_prob}")

    print("\n--- Final Answer ---")
    print(f"The final probability is {total_prob}, or {float(total_prob):.2f}.")

if __name__ == '__main__':
    solve_gossip_probability()