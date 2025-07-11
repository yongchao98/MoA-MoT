def solve():
    """
    Calculates the number of involutions in the group PSU(4, 997).
    """
    q = 997

    # An involution in PSU(4, q) is an element gZ (where Z is the center of SU(4, q))
    # such that (gZ)^2 = Z and gZ != Z. This means g^2 is in Z and g is not in Z.
    # The center Z of SU(4, 997) is {I, -I}, with size 2.
    # So we need to count elements g in SU(4, 997) where g^2 = I or g^2 = -I,
    # and g is not I or -I.

    # Case 1: g^2 = I (g is a non-central involution in SU(4, q))
    # The number of such elements is given by the formula: q^4 * (q^2 - q + 1) * (q^2 + 1)
    s1_term1 = q**4
    s1_term2 = q**2 - q + 1
    s1_term3 = q**2 + 1
    s1 = s1_term1 * s1_term2 * s1_term3

    # Case 2: g^2 = -I
    # For q = 997 = 1 (mod 4), the number of such elements is given by the formula:
    # (q^5 * (q^3 + 1) * (q^4 - 1)) / 2
    s2_term1 = q**5
    s2_term2 = q**3 + 1
    s2_term3 = q**4 - 1
    # Use integer division since the result must be an integer
    s2 = (s2_term1 * s2_term2 * s2_term3) // 2

    # The total number of pre-images in SU(4, q) is s1 + s2.
    total_s = s1 + s2

    # Each involution in PSU(4, q) corresponds to 2 elements in SU(4, q).
    # So, the number of involutions is (s1 + s2) / 2.
    result = total_s // 2

    # Print the final equation with all the numbers as requested.
    print("The number of involutions in PSU(4, 997) is calculated as follows:")
    print("Let q = 997.")
    print("\nStep 1: Count elements g in SU(4, q) with g^2 = I (non-central).")
    print("Formula: S1 = q^4 * (q^2 - q + 1) * (q^2 + 1)")
    print(f"S1 = {s1_term1} * {s1_term2} * {s1_term3}")
    print(f"S1 = {s1}")

    print("\nStep 2: Count elements g in SU(4, q) with g^2 = -I.")
    print("Formula: S2 = (q^5 * (q^3 + 1) * (q^4 - 1)) / 2")
    print(f"S2 = ({s2_term1} * {s2_term2} * {s2_term3}) / 2")
    print(f"S2 = {s2}")

    print("\nStep 3: Calculate the total number of involutions in PSU(4, 997).")
    print("Formula: Result = (S1 + S2) / 2")
    print(f"Result = ({s1} + {s2}) / 2")
    print(f"Result = {total_s} / 2")
    print(f"Result = {result}")

solve()