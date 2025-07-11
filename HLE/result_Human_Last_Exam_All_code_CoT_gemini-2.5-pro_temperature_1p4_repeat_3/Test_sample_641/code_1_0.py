def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The formula for the number of involutions in PSU(n, q) for n=4 and q === 1 (mod 4)
    # is (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2.

    # Calculate the components of the numerator
    q_pow_4 = q**4
    term1 = q**2 - q + 1
    term2 = q**2 + 1

    # The numerator is the number of involutions in SU(4, q)
    numerator = q_pow_4 * term1 * term2

    # The result is the numerator divided by the size of the center, which is 2.
    result = numerator // 2

    # As requested, printing the final equation with each number.
    print(f"The number of involutions is given by the equation:")
    print(f"({q_pow_4} * {term1} * {term2}) / 2 = {result}")

solve()