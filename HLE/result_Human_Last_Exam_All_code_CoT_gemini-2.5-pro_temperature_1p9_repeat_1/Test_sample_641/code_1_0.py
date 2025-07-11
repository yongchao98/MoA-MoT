def solve_involutions():
    """
    This function calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # Case 1: Involutions arising from elements A in SU(4,q) with A^2 = I.
    # The number of such involutions in PSU(4,q) is given by the formula:
    # n1 = (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2
    # Python's integer division // is used as the result is always an integer.
    num1 = (q**4 * (q**2 - q + 1) * (q**2 + 1)) // 2

    # Case 2: Involutions arising from elements A in SU(4,q) with A^2 = -I.
    # This case contributes involutions because q=997 is congruent to 1 (mod 4).
    # The number of such involutions in PSU(4,q) is given by the formula:
    # n2 = (q^5 * (q + 1)^2 * (q^2 - q + 1) * (q^2 + 1)) / 2
    num2 = (q**5 * (q + 1)**2 * (q**2 - q + 1) * (q**2 + 1)) // 2

    # The total number of involutions is the sum of these two cases.
    total_involutions = num1 + num2

    # Print the final equation with the computed numbers, as requested.
    print(f"{num1} + {num2} = {total_involutions}")

solve_involutions()