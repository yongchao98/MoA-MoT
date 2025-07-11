def solve_scheme_properties():
    """
    This function determines and prints the lexicographically ordered list of
    all maximal subsets of {A,B,C,D,E} for which a scheme X exists.

    The logic is based on known properties and contradictions in algebraic geometry.
    - A: dim 1 over C
    - B: projective variety over C
    - C: not reduced
    - D: affine scheme
    - E: not separated

    Forbidden combinations: {B,C}, {B,E}, {D,E}, {A,B,D}.

    The maximal subsets satisfying these constraints are determined by logical deduction:
    - {A, B}: e.g., the projective line P^1. It is maximal because adding C, D, or E leads to a contradiction.
    - {A, C, D}: e.g., Spec(C[x] x C[y]/(y^2)). It is maximal because adding B or E leads to a contradiction.
    - {A, C, E}: e.g., a non-reduced line with a doubled origin. It is maximal because adding B or D leads to a contradiction.
    - {B, C, D}: e.g., a fat point Spec(C[t]/(t^2)). It is maximal because it must have dim 0 (ruling out A) and cannot be non-separated (ruling out E).
    """

    # List of the maximal subsets, pre-sorted for lexicographical output.
    maximal_subsets = [
        ["A", "B"],
        ["A", "C", "D"],
        ["A", "C", "E"],
        ["B", "C", "D"]
    ]

    # Format each subset into the format {X,Y,Z}
    formatted_subsets = []
    for subset in maximal_subsets:
        # The elements within each subset are already alphabetically sorted.
        formatted_subsets.append("{" + ",".join(subset) + "}")

    # Join the formatted subsets into a single string.
    final_output = ", ".join(formatted_subsets)

    print(final_output)

solve_scheme_properties()