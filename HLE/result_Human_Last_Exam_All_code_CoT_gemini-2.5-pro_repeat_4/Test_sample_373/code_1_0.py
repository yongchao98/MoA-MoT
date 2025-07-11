import math

def solve():
    """
    Calculates the probability of sampling the bottleneck edge in a
    barbell graph with 10 nodes using randomized uniform gossiping.
    """
    # 1. Define the graph parameters
    total_nodes = 10

    # A standard barbell graph has two identical cliques.
    # The size of each clique is half the total number of nodes.
    clique_size = total_nodes / 2

    # 2. Calculate the degree of the bridge nodes
    # The bottleneck edge connects one "bridge node" from each clique.
    # The degree of a bridge node is its connections within its own clique (clique_size - 1)
    # plus its one connection across the bottleneck edge.
    # Degree = (clique_size - 1) + 1 = clique_size
    degree_of_bridge_node = clique_size

    # 3. Calculate the probability of the two scenarios
    # In randomized uniform gossiping, we first pick a node from all N nodes,
    # then that node picks one of its neighbors.

    # The probability of picking any single node is 1 / total_nodes.
    # The probability of a bridge node picking the other bridge node is 1 / its degree.

    # Scenario 1: Pick bridge node A, which then picks bridge node B
    prob_scenario_1 = (1 / total_nodes) * (1 / degree_of_bridge_node)

    # Scenario 2: Pick bridge node B, which then picks bridge node A
    prob_scenario_2 = (1 / total_nodes) * (1 / degree_of_bridge_node)

    # 4. The total probability is the sum of the two scenarios
    total_probability = prob_scenario_1 + prob_scenario_2

    print("Step-by-step calculation for the probability of sampling the bottleneck edge:")
    print(f"Total nodes (N): {int(total_nodes)}")
    print(f"Size of each clique (m): {int(clique_size)}")
    print(f"Degree of a bridge node (deg): (m-1) + 1 = {int(degree_of_bridge_node)}")
    print("\nThe probability is the sum of two scenarios (picking bridge node A first, or picking B first):")
    print("P = P(pick A) * P(A picks B) + P(pick B) * P(B picks A)")
    print("The final equation is:")
    print(f"P = (1 / {int(total_nodes)}) * (1 / {int(degree_of_bridge_node)}) + (1 / {int(total_nodes)}) * (1 / {int(degree_of_bridge_node)})")
    print(f"P = {prob_scenario_1} + {prob_scenario_2}")
    print(f"Total Probability = {total_probability}")
    print(f"As a fraction, this is 1/{int(1/total_probability)}")


solve()