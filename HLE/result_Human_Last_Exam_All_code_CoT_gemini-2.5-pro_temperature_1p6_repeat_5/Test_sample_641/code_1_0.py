def calculate_involutions_in_psu_4_q():
    """
    Calculates the number of involutions in the Projective Special Unitary group PSU(4, 997).

    The number of involutions in PSU(4, q) for q odd is the sum of the sizes of two distinct
    conjugacy classes of involutions.
    """
    q = 997

    # The size of the first conjugacy class of involutions (from elements g with g^2 = I)
    k1 = q**4 * (q**2 + 1) * (q**2 - q + 1)

    # The size of the second conjugacy class of involutions (from elements g with g^2 = -I)
    k2 = q**4 * (q**2 - 1) * (q**3 + 1)

    # The total number of involutions is the sum of the sizes of these two classes.
    total_involutions = k1 + k2

    # As requested, output each number in the final equation.
    print(f"The number of involutions in PSU(4,997) is the sum of the sizes of two classes.")
    print(f"Class 1 size: {k1}")
    print(f"Class 2 size: {k2}")
    print(f"\nTotal number of involutions is:\n{k1} + {k2} = {total_involutions}")

calculate_involutions_in_psu_4_q()