import sys

def solve_gossiping_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes using randomized uniform gossiping.
    """

    # 1. Define the parameters of the barbell graph
    total_nodes = 10
    # A 10-node barbell graph consists of two K5 cliques
    clique_size = total_nodes // 2

    # 2. Determine the degree of the nodes connected by the bottleneck edge
    # The degree of a node in a complete graph of size `m` is `m-1`.
    # The two nodes on the bottleneck edge have one additional connection.
    # So, their degree is (clique_size - 1) + 1.
    degree_of_bridge_node = (clique_size - 1) + 1

    # 3. Calculate the probability
    # The process: a node is chosen uniformly at random, then it chooses a neighbor uniformly at random.
    # Let the bottleneck edge be (u, v).
    # P(sample edge) = P(select u) * P(u chooses v) + P(select v) * P(v chooses u)

    # Probability of selecting any single node from the total
    prob_select_node = 1 / total_nodes

    # Probability of a bridge node choosing its bridge neighbor
    prob_choose_bridge_neighbor = 1 / degree_of_bridge_node

    # Total probability is the sum of the two ways the edge can be sampled
    total_prob = (prob_select_node * prob_choose_bridge_neighbor) + (prob_select_node * prob_choose_bridge_neighbor)

    # 4. Print the explanation and the final result
    print("This script calculates the probability of sampling the bottleneck edge in a 10-node barbell graph.")
    print("-" * 80)
    print(f"Graph Structure:")
    print(f" - Total nodes (N): {total_nodes}")
    print(f" - Structure: Two K{clique_size} cliques connected by one bottleneck edge.")
    print("\nDegree Calculation:")
    print(f" - The two nodes connected by the bottleneck edge each have a degree of ({clique_size} - 1) + 1 = {degree_of_bridge_node}.")
    print("\nProbability Formula:")
    print("P(sampling bottleneck) = P(select node u) * P(u chooses v) + P(select node v) * P(v chooses u)")
    print("\nCalculation:")
    # Print the equation with the actual numbers, as requested
    print(f"P(sampling bottleneck) = (1/{total_nodes}) * (1/{degree_of_bridge_node}) + (1/{total_nodes}) * (1/{degree_of_bridge_node})")
    print(f"                       = {prob_select_node:.2f} * {prob_choose_bridge_neighbor:.2f} + {prob_select_node:.2f} * {prob_choose_bridge_neighbor:.2f}")
    print(f"                       = 1/{total_nodes * degree_of_bridge_node} + 1/{total_nodes * degree_of_bridge_node}")
    print(f"                       = 2/{total_nodes * degree_of_bridge_node}")
    print(f"                       = {2 / (total_nodes * degree_of_bridge_node)}")
    print("-" * 80)
    print(f"The final probability is {total_prob}")


if __name__ == "__main__":
    solve_gossiping_probability()
