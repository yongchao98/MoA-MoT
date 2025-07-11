import math

def calculate_bottleneck_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    10-node barbell graph using randomized uniform gossiping.
    """

    # Step 1: Define the graph parameters based on the problem description.
    # A barbell graph with 10 nodes is interpreted as two K_5 cliques
    # connected by a single edge.
    num_nodes_total = 10
    nodes_per_clique = 5

    print(f"Analyzing a barbell graph with N = {num_nodes_total} nodes.")
    print(f"The graph consists of two complete graphs (K_{nodes_per_clique}) connected by one bottleneck edge.\n")

    # Step 2: Determine the degree of the nodes connected by the bottleneck edge.
    # Let the bottleneck edge be (u, v).
    # Degree of u = (connections inside its K_5) + (the bottleneck connection)
    degree_bridge_node = (nodes_per_clique - 1) + 1
    
    print("The gossiping process is as follows:")
    print("1. A node 'u' is picked uniformly at random from N nodes.")
    print("2. A neighbor of 'u' is picked uniformly at random from its deg(u) neighbors.\n")

    print("The two nodes connected by the bottleneck edge have a degree calculated as:")
    print(f"Degree = (Nodes in clique - 1) + 1 = ({nodes_per_clique} - 1) + 1 = {degree_bridge_node}\n")

    # Step 3: Calculate the probability of sampling the bottleneck edge (u,v).
    # P(sampling) = P(picking u) * P(neighbor v | u) + P(picking v) * P(neighbor u | v)
    # P(sampling) = (1/N) * (1/deg(u)) + (1/N) * (1/deg(v))
    
    # Calculate each term of the probability equation
    prob_u_then_v = (1 / num_nodes_total) * (1 / degree_bridge_node)
    prob_v_then_u = (1 / num_nodes_total) * (1 / degree_bridge_node)
    total_prob = prob_u_then_v + prob_v_then_u

    # For displaying the fraction
    numerator = 2
    denominator = num_nodes_total * degree_bridge_node
    common_divisor = math.gcd(numerator, denominator)
    simple_numerator = numerator // common_divisor
    simple_denominator = denominator // common_divisor

    print("The probability of sampling the bottleneck edge is calculated as:")
    print("P(edge) = P(pick first node) * P(pick second node) + P(pick second node) * P(pick first node)")
    print(f"P(edge) = (1 / {num_nodes_total}) * (1 / {degree_bridge_node}) + (1 / {num_nodes_total}) * (1 / {degree_bridge_node})")
    print(f"P(edge) = (1 / {denominator}) + (1 / {denominator})")
    print(f"P(edge) = {numerator} / {denominator}")
    print(f"P(edge) = {simple_numerator} / {simple_denominator}\n")
    print(f"The final probability is {total_prob}")


if __name__ == '__main__':
    calculate_bottleneck_probability()