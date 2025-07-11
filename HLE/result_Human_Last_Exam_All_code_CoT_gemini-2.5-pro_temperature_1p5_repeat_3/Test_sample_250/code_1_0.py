def solve_scheme_properties():
    """
    This function determines and prints the lexicographically ordered list of
    all maximal subsets of {A,B,C,D,E} for which a scheme X exists.

    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated
    """

    # Based on the analysis of contradictions between the properties,
    # the maximal subsets are pre-determined.
    # 1. {A, B, C}: Non-reduced projective curve. (e.g., V(x^2*y) in P^2)
    # 2. {A, C, D}: Non-reduced affine curve. (e.g., Spec(C[x,y]/(y^2)))
    # 3. {A, C, E}: Non-reduced non-separated curve. (e.g., non-reduced line with doubled origin)
    # 4. {B, C, D}: 0-dim, non-reduced, affine, and projective scheme. (e.g., Spec(C[eps]/eps^2))
    maximal_subsets = [
        ['A', 'B', 'C'],
        ['A', 'C', 'D'],
        ['A', 'C', 'E'],
        ['B', 'C', 'D']
    ]

    # Format the list for output as specified in the problem.
    # e.g., "{A,B,C}, {A,C,D}, ..."
    formatted_sets = []
    for s in maximal_subsets:
        formatted_sets.append(f"{{{','.join(s)}}}")

    # Print the final comma-separated list.
    print(', '.join(formatted_sets))

solve_scheme_properties()