def solve_psu_involutions():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # Part 1: Calculate the number of involutions in SU(4, q), denoted i(SU(4,q)).
    # An involution t in SU(4,q) has eigenvalues +-1 and det(t)=1,
    # which implies the (-1)-eigenspace has even dimension. For n=4, this dimension is 2.
    # The formula for the number of such elements is i(SU(4,q)) = q^4 * (q^2 - q + 1) * (q^2 + 1).
    i_su4q = q**4 * (q**2 - q + 1) * (q**2 + 1)

    # Part 2: Calculate the number of elements g in SU(4, q) such that g^2 = -I, denoted s_2(SU(4,q)).
    # These elements also correspond to involutions in the quotient group PSU(4, q).
    # The formula is s_2(SU(4,q)) = q^2 * (q^3 + 1).
    s2_su4q = q**2 * (q**3 + 1)

    # Part 3: Calculate the number of involutions in PSU(4, q).
    # The center Z of SU(4, q) has size d = gcd(4, q+1) = gcd(4, 998) = 2.
    # Involutions in PSU(4,q) arise from elements g in SU(4,q) where g^2 is in Z, and g is not in Z.
    # The total number of such pre-images is (i(SU(4,q)) + 1) + s_2(SU(4,q)).
    # We subtract the size of the center |Z|=2 (since these map to the identity) and divide by |Z|=2.
    # Number of involutions = (i(SU(4,q)) + s_2(SU(4,q)) + 1 - 2) / 2
    num_involutions = (i_su4q + s2_su4q - 1) // 2

    # --- Outputting the results step-by-step ---
    print(f"The number of involutions in PSU(4, 997) is calculated as follows, with q = {q}:\n")

    print("Step 1: Calculate the number of involutions in SU(4, q), i(SU(4,q)).")
    print("Formula: i(SU(4,q)) = q^4 * (q^2 - q + 1) * (q^2 + 1)")
    print(f"i(SU(4,997)) = {q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1)")
    print(f"i(SU(4,997)) = {i_su4q}\n")

    print("Step 2: Calculate the number of elements g in SU(4, q) with g^2 = -I, s_2(SU(4,q)).")
    print("Formula: s_2(SU(4,q)) = q^2 * (q^3 + 1)")
    print(f"s_2(SU(4,997)) = {q}^2 * ({q}^3 + 1)")
    print(f"s_2(SU(4,997)) = {s2_su4q}\n")

    print("Step 3: Combine results to find the number of involutions in PSU(4, q).")
    print("Formula: (i(SU(4,q)) + s_2(SU(4,q)) - 1) / 2")
    print(f"Number = ({i_su4q} + {s2_su4q} - 1) / 2")
    print(f"Number = {i_su4q + s2_su4q - 1} / 2")
    print(f"The total number of involutions in PSU(4, 997) is: {num_involutions}")

solve_psu_involutions()
print(f'<<<{ (997**4 * (997**2 - 997 + 1) * (997**2 + 1) + 997**2 * (997**3 + 1) - 1) // 2 }>>>')