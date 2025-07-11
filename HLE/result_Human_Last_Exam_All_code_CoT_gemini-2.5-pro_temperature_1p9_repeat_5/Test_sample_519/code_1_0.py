def solve():
    """
    This function calculates and prints the properties of the three given
    categories fibered in groupoids.
    """

    # Properties for X1 = Hilb_11(A^3)
    # Type: Scheme (S) - Representable by a scheme.
    # Separated (s) - Hilbert schemes are separated.
    # Not universally closed (uc) - Base A^3 is not proper.
    # Not irreducible (irr) - Hilb_d(A^n) is reducible for n>=3, d>=4.
    # Dimension (dim) = 33 - The component of distinct points has dim = n*d = 3*11 = 33.
    X1_props = ["S", "s", "33"]

    # Properties for X2 = [(A^4 \ V(xy-zw))/C*]
    # Type: Deligne-Mumford stack (DM) - Quotient by C* with finite stabilizers.
    # Separated (s) - Quotient of a separated scheme by a reductive group.
    # Not universally closed (uc) - Base space is not proper.
    # Irreducible (irr) - A^4 \ V(...) is irreducible.
    # Dimension (dim) = 3 - Calculated as dim(A^4) - dim(C*) = 4 - 1 = 3.
    X2_props = ["DM", "s", "irr", "3"]

    # Properties for X3 = Pic(C_0) for g=7 curve C_0
    # Type: Scheme (S) - Representable by the Picard scheme.
    # Separated (s) - It is a group scheme.
    # Not universally closed (uc) - Infinite disjoint union of proper components is not proper.
    # Not irreducible (irr) - Has infinitely many connected components.
    # Dimension (dim) = 7 - Dimension is the genus of the curve, g = 7.
    X3_props = ["S", "s", "7"]

    # Format the properties into the required string representation
    x1_str = f"[{','.join(X1_props)}]"
    x2_str = f"[{','.join(X2_props)}]"
    x3_str = f"[{','.join(X3_props)}]"

    # Combine the strings and print the final result
    final_answer = f"{x1_str} {x2_str} {x3_str}"
    print(final_answer)

solve()