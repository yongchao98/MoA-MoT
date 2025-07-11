def find_maximal_scheme_subsets():
    """
    This function prints the lexicographically ordered list of all maximal subsets of
    properties {A,B,C,D,E} that a scheme X can satisfy.

    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated

    Based on a careful analysis of the compatibility of these properties in algebraic
    geometry, the following four maximal subsets are identified:
    - {A, B, C}: e.g., a non-reduced projective curve like V(Z^2) in P^2.
    - {A, C, D}: e.g., a non-reduced affine curve like Spec(C[x,y]/(y^2)).
    - {A, C, E}: e.g., a non-reduced, non-separated curve like a "fat line with a doubled axis".
    - {B, C, D}: e.g., a 0-dimensional, non-reduced scheme like Spec(C[x]/(x^2)).

    This script formats and prints this result.
    """
    maximal_subsets = [
        "{A, B, C}",
        "{A, C, D}",
        "{A, C, E}",
        "{B, C, D}"
    ]

    print(", ".join(maximal_subsets))

find_maximal_scheme_subsets()