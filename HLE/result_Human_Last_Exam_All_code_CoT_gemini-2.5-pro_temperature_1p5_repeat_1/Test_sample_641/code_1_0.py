def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The formula for the number of involutions in PSU(4, q) for q=1 (mod 4) is:
    # N = (q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2
    # We calculate the terms of the equation first.

    q_pow_2 = q**2
    q_pow_4 = q**4

    term1 = q_pow_2 + 1
    term2 = q_pow_2 - q + 1

    numerator = q_pow_4 * term1 * term2
    result = numerator // 2

    # Print the equation with the computed values, as requested.
    print("The formula for the number of involutions is:")
    print(f"N = (q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2")
    print("\nFor q = 997:")
    print(f"N = ({q}^4 * ({q}^2 + 1) * ({q}^2 - {q} + 1)) / 2")
    print(f"N = ({q_pow_4} * {term1} * {term2}) / 2")
    print(f"N = {numerator} / 2")
    print(f"\nThe number of involutions is:")
    print(result)

solve()