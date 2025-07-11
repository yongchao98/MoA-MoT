def solve_scheme_properties():
    """
    This function determines and prints the lexicographically ordered list of all
    maximal subsets of properties {A,B,C,D,E} that a scheme X can satisfy.

    The properties are:
    A) X has dimension 1 over C
    B) X is a projective variety over C
    C) X is not reduced
    D) X is an affine scheme
    E) X is not separated

    The reasoning is as follows:
    1.  We identify fundamental contradictions between properties:
        - {A, B, D}: A scheme over a field is projective (B) and affine (D) iff it is 0-dimensional. This contradicts dimension 1 (A).
        - {B, E}: A projective scheme (B) is always separated, which contradicts being not separated (E).
        - {D, E}: An affine scheme (D) is always separated, which contradicts being not separated (E).

    2.  We argue that any maximal set must contain property C (not reduced).
        For any scheme X realizing a set of properties S from {A, B, D, E}, one can
        construct a new scheme X' = X x Spec(C[eps]/eps^2). This new scheme X' is
        not reduced (has property C) and preserves all the properties in S that X had.
        This implies that any set of properties S without C cannot be maximal, because
        S U {C} is also a valid set of properties.

    3.  Therefore, we search for maximal subsets of {A, B, D, E} that do not contain
        any of the contradictory sets identified in step 1. These are:
        - {A, B}
        - {A, D}
        - {A, E}
        - {B, D} (Note: {B, D} is possible for 0-dimensional schemes, so A is excluded).

    4.  Adding C to each of these gives the maximal sets for the original problem:
        - {A, B, C}: Example: A non-reduced projective curve, like the plane curve defined by (F(x,y,z))^2=0 where F is an irreducible homogeneous polynomial.
        - {A, C, D}: Example: A non-reduced affine curve, like Spec(C[x, eps]/(eps^2)).
        - {A, C, E}: Example: A non-reduced non-separated curve, like the "affine line with a doubled origin" constructed from non-reduced rings.
        - {B, C, D}: Example: A 0-dimensional, non-reduced scheme, like Spec(C[t]/(t^2)). It is affine, projective, and not reduced, but does not have dimension 1.

    5.  Finally, we order these sets lexicographically.
    """
    
    # The maximal subsets of {A,B,C,D,E} are determined by the logical analysis above.
    # The code below simply prints the final, pre-computed result in the required format.
    
    set1 = "{A,B,C}"
    set2 = "{A,C,D}"
    set3 = "{A,C,E}"
    set4 = "{B,C,D}"
    
    # The sets are already in lexicographical order.
    # We print them in a single line, separated by commas.
    print(f"{set1}, {set2}, {set3}, {set4}")

solve_scheme_properties()