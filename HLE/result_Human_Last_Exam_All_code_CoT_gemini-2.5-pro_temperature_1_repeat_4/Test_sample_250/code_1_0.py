import itertools

def solve_scheme_properties():
    """
    This function determines the maximal subsets of compatible properties for a scheme X.
    The logic is based on known incompatibilities between the properties.
    """
    # Define the set of properties
    properties = ['A', 'B', 'C', 'D', 'E']

    # Define the sets of incompatible properties based on algebraic geometry theorems.
    # A scheme cannot simultaneously have the properties in any of these sets.
    incompatibilities = [
        {'B', 'C'},  # A projective variety is reduced.
        {'B', 'E'},  # A projective scheme is separated.
        {'D', 'E'},  # An affine scheme over a field is separated.
        {'A', 'B', 'D'} # A scheme that is both projective and affine over a field has dimension 0.
    ]

    # Generate all non-empty subsets of the properties
    all_subsets = []
    for i in range(1, len(properties) + 1):
        for subset in itertools.combinations(properties, i):
            all_subsets.append(set(subset))

    # Filter for subsets that are "possible" (do not contain any incompatible set)
    possible_subsets = []
    for s in all_subsets:
        is_possible = True
        for inc in incompatibilities:
            if inc.issubset(s):
                is_possible = False
                break
        if is_possible:
            possible_subsets.append(s)

    # From the list of possible subsets, find the ones that are maximal.
    # A subset is maximal if it's not a proper subset of any other possible subset.
    maximal_subsets = []
    for s1 in possible_subsets:
        is_maximal = True
        for s2 in possible_subsets:
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # To ensure a canonical order for output, we sort the results.
    # First, sort the letters within each subset.
    sorted_maximal_subsets = [sorted(list(s)) for s in maximal_subsets]
    # Then, sort the list of subsets lexicographically.
    sorted_maximal_subsets.sort()

    # Format the output string as per the example.
    output_strings = [f"{{{', '.join(subset)}}}" for subset in sorted_maximal_subsets]
    final_output = ', '.join(output_strings)

    print(final_output)

solve_scheme_properties()