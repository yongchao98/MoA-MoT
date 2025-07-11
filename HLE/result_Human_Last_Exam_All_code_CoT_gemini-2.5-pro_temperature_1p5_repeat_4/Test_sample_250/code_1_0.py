import itertools

def find_maximal_subsets():
    """
    This function finds all maximal subsets of properties that a scheme can satisfy,
    based on a predefined set of rules from algebraic geometry.
    """
    
    # 1. Define the set of properties
    properties = frozenset({'A', 'B', 'C', 'D', 'E'})

    # 2. Define the pairs of incompatible properties.
    # A scheme cannot be simultaneously:
    # - B (Projective Variety) and C (Not Reduced)
    # - B (Projective, dim 1) and D (Affine, dim 1)
    # - B (Projective) and E (Not Separated)
    # - D (Affine) and E (Not Separated)
    incompatibilities = [
        frozenset({'B', 'C'}),
        frozenset({'B', 'D'}),
        frozenset({'B', 'E'}),
        frozenset({'D', 'E'}),
    ]

    # 3. Generate the power set (all non-empty subsets) of the properties
    all_subsets = []
    s = list(properties)
    for i in range(1, len(s) + 1):
        for subset in itertools.combinations(s, i):
            all_subsets.append(frozenset(subset))

    # 4. Filter out subsets that contain any incompatible pair
    consistent_subsets = []
    for subset in all_subsets:
        is_consistent = True
        for inc_pair in incompatibilities:
            if inc_pair.issubset(subset):
                is_consistent = False
                break
        if is_consistent:
            consistent_subsets.append(subset)

    # 5. From the consistent subsets, find the ones that are maximal.
    # A subset S1 is maximal if no other consistent subset S2 is a superset of S1.
    # Note: The existence of a scheme for each candidate is confirmed by the reasoning
    # provided in the text explanation.
    # {A, B}: P^1_C
    # {A, C, D}: Spec(C[x,y]/(y^2))
    # {A, C, E}: Non-reduced line with a doubled origin
    maximal_subsets = []
    for s1 in consistent_subsets:
        is_maximal = True
        for s2 in consistent_subsets:
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # 6. Sort and format the results for the final output.
    # First, convert each frozenset to a sorted list of characters for consistent ordering.
    sorted_subsets_as_lists = [sorted(list(s)) for s in maximal_subsets]
    
    # Second, sort the list of lists lexicographically.
    sorted_subsets_as_lists.sort()
    
    # Format each subset as a string like "{A, B}"
    output_strings = [f"{{{', '.join(s)}}}" for s in sorted_subsets_as_lists]

    # Join the formatted strings into the final result string
    final_output = ', '.join(output_strings)

    print(final_output)

find_maximal_subsets()