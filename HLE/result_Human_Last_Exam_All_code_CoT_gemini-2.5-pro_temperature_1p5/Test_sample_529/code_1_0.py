import collections
from itertools import chain, combinations

def solve_count_ans(graph, k):
    """
    Calculates the number of answers for the formula phi_k for a given graph G and integer k.
    An answer is a tuple of k vertices (x_1, ..., x_k) for which there exists
    a common neighbor y.

    The method uses the Principle of Inclusion-Exclusion on the distinct neighborhoods
    of the graph's vertices.
    """
    print(f"Graph vertices: {list(graph.keys())}")
    print(f"Graph edges: {[(u, v) for u, neighbors in graph.items() for v in neighbors if u < v]}")
    print(f"k = {k}\n")

    # Step 1: Find all distinct neighborhoods
    # A neighborhood N(u) is represented as a frozenset of its neighbors.
    # We map each distinct neighborhood to a list of vertices that have it.
    distinct_neighborhoods = collections.defaultdict(list)
    for vertex in sorted(graph.keys()):
        # frozenset is used because sets are mutable and cannot be dict keys
        neighborhood = frozenset(graph[vertex])
        distinct_neighborhoods[neighborhood].append(vertex)

    print("Found distinct neighborhoods:")
    named_neighborhoods = {}
    for i, (neighborhood, vertices) in enumerate(distinct_neighborhoods.items()):
        name = f"N{i+1}"
        print(f"  {name}: {set(neighborhood)} (for vertices {vertices})")
        named_neighborhoods[name] = set(neighborhood)
    print("-" * 20)
    
    # The set of answers is the union of N^k for all distinct neighborhoods N.
    # We use PIE to compute the size of this union.
    # Formula: |U N_i^k| = sum_{I non-empty} (-1)^{|I|-1} |intersect(N_i for i in I)|^k
    
    total_count = 0
    final_equation_terms = []

    distinct_list = list(named_neighborhoods.values())
    
    # Step 2: Iterate through all non-empty subsets of the set of distinct neighborhoods
    indices = range(len(distinct_list))
    for i in range(1, len(indices) + 1):
        # i is the size of the subsets of neighborhoods to intersect
        term_sum_for_size_i = 0
        
        subset_terms_strings = []

        for subset_indices in combinations(indices, i):
            # Calculate intersection of neighborhoods in the current subset
            if not subset_indices:
                continue

            # Start with the first neighborhood in the subset
            intersection = set(distinct_list[subset_indices[0]])
            # Intersect with the rest
            for j in range(1, len(subset_indices)):
                intersection.intersection_update(distinct_list[subset_indices[j]])
            
            # Calculate the term for the PIE formula
            term = len(intersection) ** k
            if i % 2 == 0: # if size of subset is even, subtract
                term_sum_for_size_i -= term
                subset_terms_strings.append(f"{len(intersection)}^{k}")
            else: # if size of subset is odd, add
                term_sum_for_size_i += term
                subset_terms_strings.append(f"{len(intersection)}^{k}")
        
        if subset_terms_strings:
            op = "+" if i % 2 != 0 else "-"
            # Join all |N_i cap N_j ...|^k terms for this subset size
            sum_str = " + ".join(subset_terms_strings)
            # Add brackets for sums of more than one term
            if len(subset_terms_strings) > 1:
                final_equation_terms.append(f"{op} ({sum_str})")
            else:
                final_equation_terms.append(f"{op} {sum_str}")

        total_count += term_sum_for_size_i
    
    print("Calculation using Principle of Inclusion-Exclusion:")
    # The first term does not need a leading '+'
    if final_equation_terms:
        first_term = final_equation_terms[0]
        if first_term.startswith("+ "):
            final_equation_terms[0] = first_term[2:]

    final_equation_str = " ".join(final_equation_terms)
    print(f"Count = {final_equation_str}")
    print(f"      = {total_count}")

# Example Usage:
# Adjacency list representation of a graph. Using a path graph on 4 vertices.
# 1 -- 2 -- 3 -- 4
example_graph = {
    1: {2},
    2: {1, 3},
    3: {2, 4},
    4: {3}
}
example_k = 2

solve_count_ans(example_graph, example_k)