import math

def solve_clustering_coefficient():
    """
    This function calculates the clustering coefficient of the central node
    of the recursively defined graph.
    """

    # Step 1 & 2: Define recurrence relations for nodes N(L) and edges M(L).
    # Level L=0 (the initial blue K4)
    # A K4 has 4 nodes and 6 edges.
    n0 = 4
    m0 = 6

    # Level L=1 (the 16-node structure)
    # N(L+1) = 4 * N(L)
    # M(L+1) = 4 * M(L) + 6 * N(L)
    n1 = 4 * n0
    m1 = 4 * m0 + 6 * n0

    # Step 3: Calculate properties of the final Level L=2 graph (64 nodes).
    n2 = 4 * n1
    m2 = 4 * m1 + 6 * n1

    print(f"The final graph is a Level 2 structure.")
    print(f"It has a total of N = {n2} nodes and M = {m2} edges.")
    print("-" * 20)

    # Step 4: Determine k and E for the central node.
    # The central node is connected to all other nodes in the graph.
    # Its degree 'k' is the total number of nodes minus one.
    k = n2 - 1

    # The number of edges 'E' between its neighbors is the total number of edges
    # in the graph minus the edges connected to the central node (which is its degree k).
    E = m2 - k

    print(f"The central node is connected to all other nodes.")
    print(f"Its degree (number of neighbors) is k = N - 1 = {n2} - 1 = {k}.")
    print(f"The number of edges between its neighbors is E = M - k = {m2} - {k} = {E}.")
    print("-" * 20)

    # Step 5: Calculate the clustering coefficient C.
    # C = 2 * E / (k * (k - 1))
    numerator = 2 * E
    denominator = k * (k - 1)
    
    # Simplify the fraction for a cleaner representation
    common_divisor = math.gcd(numerator, denominator)
    simplified_numerator = numerator // common_divisor
    simplified_denominator = denominator // common_divisor

    # Calculate the decimal value
    c = numerator / denominator

    print("The clustering coefficient C is calculated as:")
    print(f"C = (2 * E) / (k * (k - 1))")
    print(f"C = (2 * {E}) / ({k} * ({k} - 1))")
    print(f"C = {numerator} / {denominator}")
    print(f"As a simplified fraction, C = {simplified_numerator}/{simplified_denominator}")
    print(f"As a decimal, C ≈ {c}")

if __name__ == "__main__":
    solve_clustering_coefficient()