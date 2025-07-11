import math

def calculate_clustering_coefficient():
    """
    Calculates the clustering coefficient of the central node of the 64-node graph.
    """
    # Step 1: Define the degree (k) of the central node.
    # The graph has 64 nodes, and the central node connects to all others.
    k = 64 - 1

    # Step 2: Define the number of edges between the neighbors (E).
    # As derived in the explanation, this is 174.
    E = 174

    # Step 3: Calculate the clustering coefficient C.
    # C = 2 * E / (k * (k - 1))
    numerator = 2 * E
    denominator = k * (k - 1)
    
    C = numerator / denominator

    # Simplify the fraction for a cleaner representation
    common_divisor = math.gcd(numerator, denominator)
    simple_numerator = numerator // common_divisor
    simple_denominator = denominator // common_divisor

    print("The formula for the clustering coefficient is C = 2 * E / (k * (k - 1))")
    print(f"For the central node, the degree k = {k}")
    print(f"The number of edges between its neighbors is E = {E}")
    print("\nPlugging in the numbers:")
    print(f"C = (2 * {E}) / ({k} * ({k} - 1))")
    print(f"C = {numerator} / {denominator}")
    print(f"As a simplified fraction, C = {simple_numerator} / {simple_denominator}")
    print(f"As a decimal, C is approximately: {C}")

calculate_clustering_coefficient()