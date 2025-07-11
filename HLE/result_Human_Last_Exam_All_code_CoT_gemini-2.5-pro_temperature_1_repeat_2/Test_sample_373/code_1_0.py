import math

def calculate_bottleneck_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a barbell graph
    with 10 nodes using randomized uniform gossiping.
    """
    # 1. Define the Barbell Graph parameters
    total_nodes = 10

    # A barbell graph with 2m nodes consists of two K_m cliques.
    # For 10 nodes, 2*m = 10, so m = 5.
    # The graph consists of two K_5 cliques.
    m = total_nodes // 2

    # 2. Determine the degree of the nodes connected by the bottleneck edge.
    # A node in a K_m clique is connected to (m-1) other nodes in the clique.
    # The bottleneck edge adds 1 more connection to these specific nodes.
    # So, degree = (m - 1) + 1 = m.
    degree_bottleneck_node = m

    # 3. Calculate the probability using the uniform gossiping formula.
    # P(u, v) = (1/N) * (1/deg(u)) + (1/N) * (1/deg(v))
    # where N is the total number of nodes, and u and v are the bottleneck nodes.
    prob = (1 / total_nodes) * (1 / degree_bottleneck_node) + \
           (1 / total_nodes) * (1 / degree_bottleneck_node)

    # 4. Print the steps and the final equation
    print("A barbell graph with 10 nodes has two complete graphs (K_5) connected by a bottleneck edge.")
    print(f"The total number of nodes (N) is {total_nodes}.")
    print(f"The degree of each node connected by the bottleneck edge is {degree_bottleneck_node}.")
    print("\nThe probability P of sampling the bottleneck edge (u,v) is calculated as:")
    print("P = (1/N) * (1/degree(u)) + (1/N) * (1/degree(v))")
    print(f"P = (1/{total_nodes}) * (1/{degree_bottleneck_node}) + (1/{total_nodes}) * (1/{degree_bottleneck_node})")
    print(f"P = 1/{total_nodes * degree_bottleneck_node} + 1/{total_nodes * degree_bottleneck_node}")
    print(f"P = 2/{total_nodes * degree_bottleneck_node}")
    print(f"P = {prob}")
    print(f"As a fraction, the probability is 1/{int(1/prob)}.")


calculate_bottleneck_probability()
<<<0.04>>>