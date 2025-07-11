def describe_filled_groups():
    """
    This function describes the nonabelian filled groups of order 2*q^m,
    where q is an odd prime and m is a natural number. The classification
    is based on established results in group theory.
    """

    print("The nonabelian filled groups G of order 2*q^m (q is an odd prime, m >= 1) are classified as follows.")
    print("They are all semidirect products of their Sylow q-subgroup P and a cyclic group of order 2, C_2.")
    print("The classification depends on the structure of P.\n")

    # Case 1: P is cyclic
    print("Case 1: The Sylow q-subgroup P is cyclic, P = C_{q^m}")
    print("=====================================================")
    print("Description: For any odd prime q and integer m >= 1, the group is the dihedral group D_{2*q^m}.")
    print("Structure: G = C_{q^m} : C_2")
    print("Presentation: <r, s | r^(q^m) = 1, s^2 = 1, srs = r^(-1)>")
    print("Example (q=3, m=1): The group is D_6 of order 6. This group is also the symmetric group S_3.")
    print("Example (q=5, m=2): The group is D_50 of order 50.\n")

    # Case 2: P is elementary abelian
    print("Case 2: The Sylow q-subgroup P is elementary abelian, P = (C_q)^m")
    print("==================================================================")
    print("Description: For any odd prime q and integer m >= 1, the group is a semidirect product G = (C_q)^m : C_2 where the C_2 acts by inversion on P (every element is mapped to its inverse).")
    print("Structure: G = (C_q)^m : C_2 where b*x*b^(-1) = x^(-1) for x in P, b in C_2.")
    print("Example (q=5, m=2): The group is (C_5 x C_5) : C_2 of order 50.")
    print("Note: When m=1, this group is D_{2q}, which is also covered by Case 1.\n")

    print("Exceptional Case for P = (C_q)^m:")
    print("-----------------------------------")
    print("There is one special case for q=3 and m=2 (order 18).")
    print("Group: C_3 x S_3")
    print("Structure: (C_3 x C_3) : C_2. Here the C_2 action is not inversion on the whole Sylow 3-subgroup.")
    print("The generator of C_2 acts trivially on one C_3 factor and by inversion on the other.\n")

    # Case 3: P is non-abelian
    print("Case 3: The Sylow q-subgroup P is non-abelian")
    print("=================================================")
    print("Description: This occurs only for q=3 and m=3 (order 54).")
    print("Group: P : C_2, where P is the extra-special group of order 27 and exponent 3.")
    print("P is also known as the Heisenberg group over F_3, H(F_3).")
    print("Structure: The C_2 subgroup acts by inversion on P, mapping every element of P to its inverse.")


if __name__ == "__main__":
    describe_filled_groups()