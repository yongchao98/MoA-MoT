def solve():
    """
    This function calculates the properties of three categories fibered in groupoids and prints the result.
    """
    
    # Properties for X1: Hilb_11(A^3)
    # Type: Scheme (S)
    # Separated: yes (s)
    # Universally closed: no
    # Irreducible: no
    # Dimension: 3 * 11 = 33
    x1_profile = "[S,s,33]"

    # Properties for X2: [ (A^4 \ V(xy-zw)) / C* ]
    # Type: Deligne-Mumford stack (DM) since stabilizers are finite
    # Separated: yes (s)
    # Universally closed: no
    # Irreducible: yes (irr) as it's a quotient of an irreducible space
    # Dimension: dim(A^4) - dim(C*) = 4 - 1 = 3
    x2_profile = "[DM,s,irr,3]"

    # Properties for X3: Pic(C_0) for g=7 curve C_0
    # Type: Algebraic stack (A) since Aut(L) = G_m is not finite
    # Separated: yes (s)
    # Universally closed: no, Pic(C_0) is an infinite disjoint union
    # Irreducible: no
    # Dimension: equal to the genus g = 7
    x3_profile = "[A,s,7]"
    
    final_answer = f"{x1_profile} {x2_profile} {x3_profile}"
    print(final_answer)

solve()
print("<<<[S,s,33] [DM,s,irr,3] [A,s,7]>>>")