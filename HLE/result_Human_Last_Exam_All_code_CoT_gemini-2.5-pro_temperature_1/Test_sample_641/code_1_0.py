def solve_psu_involutions():
    """
    Calculates the number of involutions in PSU(4, q) for a given q.
    """
    q = 997

    # The number of involutions in PSU(4, q) is given by the formula:
    # (q^4 * (q^2 - q + 1) * (q^2 + q + 2)) / 2
    # We will calculate each term of the numerator first.

    # Calculate q^4
    term1 = q**4

    # Calculate q^2 - q + 1
    term2 = q**2 - q + 1

    # Calculate q^2 + q + 2
    term3 = q**2 + q + 2

    # The full numerator is term1 * term2 * term3
    numerator = term1 * term2 * term3

    # The final result is the numerator divided by 2
    result = numerator // 2

    print("To find the number of involutions in PSU(4, 997), we use the formula:")
    print("(q^4 * (q^2 - q + 1) * (q^2 + q + 2)) / 2\n")
    print(f"For q = {q}, the terms are:")
    print(f"q^4 = {term1}")
    print(f"q^2 - q + 1 = {term2}")
    print(f"q^2 + q + 2 = {term3}\n")
    print("The final equation is:")
    print(f"({term1} * {term2} * {term3}) / 2\n")
    print(f"The number of involutions is: {result}")

solve_psu_involutions()