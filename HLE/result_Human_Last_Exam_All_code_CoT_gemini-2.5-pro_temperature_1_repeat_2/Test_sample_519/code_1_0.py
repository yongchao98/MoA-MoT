def solve():
    """
    This function provides the properties of the three given categories fibered in groupoids.
    """
    # Properties for X_1: [Type, separated, (not uc), (not irr), dimension]
    # Type=Scheme(S), separated(s), dim=33
    x1_profile = "[S, s, 33]"

    # Properties for X_2: [Type, separated, (not uc), irreducible, dimension]
    # Type=Algebraic Space(S), separated(s), irreducible(irr), dim=3
    x2_profile = "[S, s, irr, 3]"

    # Properties for X_3: [Type, separated, (not uc), (not irr), dimension]
    # Type=Algebraic Stack(A), separated(s), dim=8
    x3_profile = "[A, s, 8]"

    # Combine the profiles into a single string as requested
    final_answer = " ".join([x1_profile, x2_profile, x3_profile])
    print(final_answer)

solve()