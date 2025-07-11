def calculate_bottleneck_probability():
    """
    Calculates the probability of sampling the bottleneck edge in a
    10-node barbell graph during randomized uniform gossiping.
    """
    # Total number of nodes in the graph
    N = 10

    # Number of nodes in each complete graph (clique)
    # 2 * m = N, so m = 5
    m = N // 2

    # The degree of a bridgehead node is (m-1) connections within its
    # clique plus 1 connection for the bottleneck edge.
    degree_bridgehead = (m - 1) + 1

    # Probability of selecting a specific node (e.g., a bridgehead)
    prob_select_node = 1 / N

    # Probability of a selected bridgehead node choosing the other bridgehead node
    prob_bridgehead_chooses_neighbor = 1 / degree_bridgehead

    # The probability can be calculated from two mutually exclusive events:
    # 1. Bridgehead A is picked and it chooses bridgehead B.
    # 2. Bridgehead B is picked and it chooses bridgehead A.
    prob_one_way = prob_select_node * prob_bridgehead_chooses_neighbor
    total_prob = 2 * prob_one_way

    print("Step-by-step calculation for the probability of sampling the bottleneck edge:")
    print(f"1. Total nodes in the graph (N): {N}")
    print(f"2. Nodes in each clique (m): {m}")
    print(f"3. Degree of a bridgehead node (d_bridge = m-1 + 1): {degree_bridgehead}")
    print(f"4. Probability of selecting a bridgehead node (1/N): 1/{N}")
    print(f"5. Probability of a bridgehead node choosing the other across the bottleneck (1/d_bridge): 1/{degree_bridgehead}")
    
    print("\nThe total probability is the sum of the probabilities of two events:")
    print(" (P(select bridgehead A) * P(A chooses B)) + (P(select bridgehead B) * P(B chooses A))")
    
    # Using f-string to format the final equation string
    final_equation = (
        f"(1/{N}) * (1/{degree_bridgehead}) + (1/{N}) * (1/{degree_bridgehead}) = "
        f"1/{N*degree_bridgehead} + 1/{N*degree_bridgehead} = "
        f"2/{N*degree_bridgehead} = 1/{int((N*degree_bridgehead)/2)} = {total_prob}"
    )

    print("\nFinal Equation:")
    print(final_equation)

calculate_bottleneck_probability()
<<<0.04>>>