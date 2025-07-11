import math

def solve_gossip_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    10-node barbell graph during a uniform gossiping step.
    """
    # 1. Define the graph structure.
    # A 10-node barbell graph consists of two K_5 (5-node complete graphs)
    # connected by a single bottleneck edge.
    total_nodes = 10
    nodes_per_clique = 5

    # 2. Determine the degree of the nodes (u, v) connected by the bottleneck.
    # Degree within the clique is (nodes_per_clique - 1).
    # The bottleneck edge adds 1 to the degree.
    degree_bottleneck_node = (nodes_per_clique - 1) + 1

    # Assign variables for the formula
    N = total_nodes
    deg_u = degree_bottleneck_node
    deg_v = degree_bottleneck_node

    # 3. Calculate the probability using the uniform gossiping formula.
    # P(u,v) = P(pick u)*P(pick v|u) + P(pick v)*P(pick u|v)
    # P(u,v) = (1/N)*(1/deg(u)) + (1/N)*(1/deg(v))
    total_prob_float = (1/N) * (1/deg_u) + (1/N) * (1/deg_v)

    # Calculate the fraction for a clear explanation
    numerator = 2
    denominator = N * deg_u
    common_divisor = math.gcd(numerator, denominator)
    simple_num = numerator // common_divisor
    simple_den = denominator // common_divisor

    # 4. Print the final equation with all the numbers.
    print("The probability of sampling the bottleneck edge (u,v) is calculated as:")
    print("P(u,v) = (1/Total_Nodes) * (1/deg(u)) + (1/Total_Nodes) * (1/deg(v))")
    print("\nSubstituting the values for the 10-node barbell graph:")
    print(f"P(u,v) = (1/{N}) * (1/{deg_u}) + (1/{N}) * (1/{deg_v}) = {numerator}/{denominator} = {simple_num}/{simple_den} = {total_prob_float}")

solve_gossip_probability()