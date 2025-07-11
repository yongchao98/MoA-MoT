import collections

def find_maximal_subsets():
    """
    This function finds all maximal subsets of properties {A,B,C,D,E}
    that a scheme can simultaneously satisfy.
    """
    properties = sorted(['A', 'B', 'C', 'D', 'E'])

    # Based on algebraic geometry principles, these combinations are impossible.
    inconsistent_sets = [
        {'A', 'B', 'D'},  # Projective + Affine + Dim 1 -> Dim 0, contradiction.
        {'B', 'E'},       # Projective implies separated.
        {'D', 'E'},       # Affine implies separated.
    ]

    # Generate the power set of all properties, from largest to smallest.
    from itertools import chain, combinations
    all_subsets = list(chain.from_iterable(combinations(properties, r) for r in range(len(properties) + 1, -1, -1)))
    all_subsets = [set(s) for s in all_subsets]

    # Filter for consistent subsets
    consistent_subsets = []
    for s in all_subsets:
        is_consistent = True
        for inconsistent in inconsistent_sets:
            if inconsistent.issubset(s):
                is_consistent = False
                break
        if is_consistent:
            consistent_subsets.append(s)

    # Filter for maximal subsets from the consistent ones.
    # A consistent set s1 is maximal if no other consistent set s2 is a proper superset of s1.
    maximal_subsets = []
    for s1 in consistent_subsets:
        is_maximal = True
        for s2 in consistent_subsets:
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # Sort the final list of sets lexicographically for the output.
    # First, sort the elements within each set.
    sorted_sets_internal = [sorted(list(s)) for s in maximal_subsets]
    # Then, sort the list of sets.
    sorted_maximal_subsets = sorted(sorted_sets_internal)

    # Format the output string as per the example.
    def format_set_string(s_list):
        return '{' + ','.join(s_list) + '}'

    result_string = ', '.join([format_set_string(s) for s in sorted_maximal_subsets])
    print(result_string)

find_maximal_subsets()