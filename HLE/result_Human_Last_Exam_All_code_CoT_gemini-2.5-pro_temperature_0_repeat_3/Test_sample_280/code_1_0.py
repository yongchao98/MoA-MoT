import math

def combinations(n, k):
    """Calculates the number of combinations of k items from a set of n."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_researcher_problem():
    """
    Solves the AIME 2021 problem about machine learning researchers.
    """
    # Step 1 & 2: Define graph parameters
    N = 42  # Number of researchers (vertices)
    d = 24  # Degree of each vertex in the co-author graph G

    # Step 3: Calculate total number of possible tables (triples of vertices)
    total_triples = combinations(N, 3)

    # Step 4: Use combinatorial counting to find the sum of c_0 and c_3
    # c_0 is the number of triples with 0 edges (independent sets)
    # c_3 is the number of triples with 3 edges (triangles)

    # Total number of edges in G
    num_edges_G = N * d // 2

    # Total number of paths of length 2 (cherries) in G
    # This is the sum over all vertices of C(deg(v), 2)
    num_cherries_G = N * combinations(d, 2)

    # The sum of edges over all possible triples is |E| * (N-2)
    # This sum also equals c_1 + 2*c_2 + 3*c_3
    sum_edges_in_triples = num_edges_G * (N - 2)

    # The sum of cherries over all possible triples is the total number of cherries
    # This sum also equals c_2 + 3*c_3
    sum_cherries_in_triples = num_cherries_G

    # From the above, we can derive c_0 + c_3 without finding c_1 and c_2
    # c_0 + c_3 = total_triples - (c_1 + c_2)
    # After substitution, this simplifies to:
    c0_plus_c3 = total_triples - sum_edges_in_triples + sum_cherries_in_triples

    # Step 5: Use the given information and the symmetry argument
    # N_0 is the number of constellations with no co-authors (K3-factors of G')
    N_0 = 2027

    # We assume a symmetric relationship where N_0 = c_3 and N_3 = c_0.
    # From N_0 = c_3, we get:
    c_3 = N_0
    # From c_0 + c_3 = constant, we get:
    c_0 = c0_plus_c3 - c_3
    # The answer N_3 is assumed to be c_0
    N_3 = c_0

    # Print the final calculation step-by-step
    print("The key relationship derived from graph properties is:")
    print(f"c_0 + c_3 = C({N}, 3) - (|E| * ({N}-2)) + (N * C({d}, 2))")
    print(f"c_0 + c_3 = {total_triples} - ({num_edges_G} * {N-2}) + {num_cherries_G}")
    print(f"c_0 + c_3 = {total_triples} - {sum_edges_in_triples} + {sum_cherries_in_triples}")
    print(f"c_0 + c_3 = {c0_plus_c3}\n")

    print("Given N_0 = 2027 (constellations with no co-authors).")
    print("By symmetry, we assume N_0 = c_3 and the answer N_3 = c_0.")
    print(f"So, the number of constellations with all co-authors is:")
    print(f"N_3 = (c_0 + c_3) - c_3 = (c_0 + c_3) - N_0")
    print(f"N_3 = {c0_plus_c3} - {N_0}")
    print(f"N_3 = {N_3}")

solve_researcher_problem()