def solve():
    """
    Calculates the number of involutions in PSU(4, 997).

    The number of involutions in PSU(n, q) for n=4 and q = 1 (mod 4) is given by
    the formula: (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2.
    """
    q = 997

    # Calculate the components of the formula
    term1 = q**4
    term2 = q**2 - q + 1
    term3 = q**2 + 1

    # Calculate the numerator of the formula
    numerator = term1 * term2 * term3

    # The denominator is the size of the center of SU(4,q)
    denominator = 2

    # The result is the number of involutions
    result = numerator // denominator

    print(f"To find the number of involutions in PSU(4, {q}), we use the formula:")
    print(f"(q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2\n")
    print(f"For q = {q}:")
    print(f"The first term, q^4, is: {term1}")
    print(f"The second term, q^2 - q + 1, is: {term2}")
    print(f"The third term, q^2 + 1, is: {term3}")
    print(f"\nThe full numerator is the product of these terms:")
    print(f"Numerator = {term1} * {term2} * {term3} = {numerator}")
    print(f"The denominator is: {denominator}")
    print(f"\nThe final result is:")
    print(f"Number of involutions = {numerator} // {denominator} = {result}")

solve()