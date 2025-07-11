import itertools

def solve_count_ans(G, k):
    """
    Calculates the number of answers for the formula phi_k for a given graph G and integer k.
    An answer is a k-tuple of vertices (x_1, ..., x_k) for which there exists a common neighbor y.
    
    The method uses the Principle of Inclusion-Exclusion on the sets of k-tuples formed
    from distinct vertex neighborhoods.
    """

    # Step 1: Find distinct neighborhoods
    neighborhoods = [set(G.get(i, [])) for i in range(len(G))]
    distinct_neighborhoods_map = {}
    for nbh in neighborhoods:
        key = frozenset(nbh)
        if key not in distinct_neighborhoods_map:
            distinct_neighborhoods_map[key] = key
    
    distinct_sets = list(distinct_neighborhoods_map.values())
    L = len(distinct_sets)
    num_vertices = len(G)

    print(f"The graph has {num_vertices} vertices and {L} distinct neighborhoods.")
    for i, nbh_set in enumerate(distinct_sets):
        print(f"  N_{i+1} = {sorted(list(nbh_set))}")
    print(f"k = {k}\n")
    print("The number of answers is |Union_{i=1 to L} (N_i)^k|.")
    print("Using the Principle of Inclusion-Exclusion, the formula is:")
    print("  Sum_{I subset {1..L}, I not empty} (-1)^(|I|-1) * |Intersection_{i in I} N_i|^k\n")
    print("Calculation:")
    
    total_count = 0
    
    equation_str_parts = []
    values_str_parts = []

    for i in range(1, L + 1):
        sign = (-1)**(i - 1)
        term_sum = 0
        term_sum_str = []
        
        for subset_indices in itertools.combinations(range(L), i):
            intersection_set = set(distinct_sets[subset_indices[0]])
            for j in range(1, len(subset_indices)):
                intersection_set.intersection_update(distinct_sets[subset_indices[j]])
            
            term = len(intersection_set)**k
            term_sum += term
            term_sum_str.append(f"{len(intersection_set)}^{k}")

        total_count += sign * term_sum

        op_char = "+" if sign > 0 else "-"
        # Don't print + for the first term
        if i == 1:
            op_char = ""

        if len(term_sum_str) > 1:
            equation_str_parts.append(f"{op_char} ({' + '.join(term_sum_str)})")
        else:
            equation_str_parts.append(f"{op_char} {term_sum_str[0]}")
            
        values_str_parts.append(f"{op_char} {term_sum}")

    print("Count =", ' '.join(equation_str_parts).lstrip())
    print("      =", ' '.join(values_str_parts).lstrip())
    print(f"\nFinal Answer: {total_count}")


# Example Graph: Path P_4
# Vertices: 0, 1, 2, 3
# Edges: (0,1), (1,2), (2,3)
G_p4 = {
    0: [1],
    1: [0, 2],
    2: [1, 3],
    3: [2]
}
k_p4 = 2

solve_count_ans(G_p4, k_p4)
