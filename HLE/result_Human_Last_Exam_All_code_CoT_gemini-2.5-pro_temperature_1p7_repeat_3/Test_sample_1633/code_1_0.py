import math

def analyze_degree_constraints(n):
    """
    Analyzes the constraints on the graph's degree distribution
    to disprove statement D from the problem.

    Statement D: At least n/4 vertices must reach degree ⌈log(n)⌉
    """

    print(f"Analyzing network constraints for n = {n} nodes.")

    # In a Watts-Strogatz model with initial degree k_0 = 6,
    # the total number of edges is n * k_0 / 2.
    total_edges = n * 6 / 2
    # The sum of degrees of all vertices in any graph is 2 * |E|.
    actual_total_degree = 2 * total_edges

    print(f"The total degree of all vertices in the graph must be exactly 2 * (3*n) = {int(actual_total_degree)}.")
    print("-" * 30)

    # Now, let's calculate the minimum possible total degree under the assumption of statement D.
    print("Testing the hypothesis from Statement D...")

    num_high_degree_nodes = n / 4
    # The degree cap is used as the degree for these nodes.
    high_degree = math.ceil(math.log(n))
    degree_sum_from_high_nodes = num_high_degree_nodes * high_degree

    num_remaining_nodes = n - num_high_degree_nodes
    # The minimum degree is k_0 / 2.
    min_degree_for_any_node = 6 / 2
    min_degree_sum_from_remaining_nodes = num_remaining_nodes * min_degree_for_any_node

    # The minimum total degree if statement D were true.
    min_total_degree_if_d = degree_sum_from_high_nodes + min_degree_sum_from_remaining_nodes

    print(f"If statement D is true:")
    print(f"  Contribution from high-degree nodes:")
    print(f"  Equation: (n/4) * ceil(log(n)) = {num_high_degree_nodes} * {high_degree} = {degree_sum_from_high_nodes}")
    print(f"  Contribution from remaining nodes (at minimum possible degree):")
    print(f"  Equation: (3n/4) * 3 = {num_remaining_nodes} * {min_degree_for_any_node} = {min_degree_sum_from_remaining_nodes}")
    print("\n  The minimum possible total degree would be:")
    print(f"  Final Equation: {degree_sum_from_high_nodes} + {min_degree_sum_from_remaining_nodes} = {min_total_degree_if_d}")

    print("-" * 30)
    if min_total_degree_if_d > actual_total_degree:
        print(f"Conclusion: Statement D is FALSE.")
        print(f"Reason: The required minimum total degree ({min_total_degree_if_d}) exceeds the graph's actual total degree ({int(actual_total_degree)}).")
    else:
        # This case won't be hit for sufficiently large n
        print("Conclusion: Statement D cannot be disproven by this method for this n.")

# Use a sufficiently large n to show the asymptotic behavior.
# For log(n) > 15, (n/4)*log(n) > (15/4)*n > 3.75n. And (3n/4)*3 = 2.25n. Sum > 6n.
# log(n) > 15 means n > e^15 ~= 3.2 million. Let's use 4,000,000.
analyze_degree_constraints(4000000)
