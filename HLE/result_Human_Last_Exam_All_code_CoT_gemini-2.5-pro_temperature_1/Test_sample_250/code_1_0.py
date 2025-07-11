def find_maximal_subsets():
    """
    This function determines and prints the lexicographically ordered list of
    all maximal subsets of {A,B,C,D,E} that can be satisfied by a scheme X.

    The properties are:
    A) dim(X) = 1 over C
    B) X is a projective scheme over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated

    The impossible combinations are:
    - {B, E}: Projective schemes are separated.
    - {D, E}: Affine schemes are separated.
    - {A, B, D}: A scheme that is both projective and affine over a field has dimension 0.

    Based on this, the maximal subsets are found to be:
    - {A, B, C}: e.g., a non-reduced line in projective space.
    - {A, C, D}: e.g., a non-reduced affine line.
    - {A, C, E}: e.g., a non-reduced, non-separated line (like a line with a doubled origin, but non-reduced).
    - {B, C, D}: e.g., a non-reduced point (dim 0), which is both affine and projective.
    """
    
    # The list of maximal sets, with properties already sorted internally.
    maximal_sets = [
        ['A', 'B', 'C'],
        ['A', 'C', 'D'],
        ['A', 'C', 'E'],
        ['B', 'C', 'D']
    ]

    # The list is already lexicographically sorted.
    # Now, format it for printing as per the example.
    
    output_parts = []
    for s in maximal_sets:
        # Format each set as "{X,Y,Z}"
        output_parts.append("{" + ",".join(s) + "}")
    
    # Join the formatted sets with ", "
    final_output = ", ".join(output_parts)
    
    print(final_output)

find_maximal_subsets()