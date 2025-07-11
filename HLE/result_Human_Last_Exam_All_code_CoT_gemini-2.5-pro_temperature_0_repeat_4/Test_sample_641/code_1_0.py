def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The number of involutions in PSU(4, q) for q=1 (mod 4) is
    # (1/2) * q^4 * (q^2 + 1) * (q^2 - q + 1)

    q_2 = q * q
    q_4 = q_2 * q_2

    term1 = q_2 + 1
    term2 = q_2 - q + 1

    # Use integer division
    result = (q_4 * term1 * term2) // 2

    print(f"The number of involutions is given by the formula:")
    print(f"(q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2")
    print(f"Substituting q = {q}:")
    print(f"({q_4} * {term1} * {term2}) // 2 = {result}")

solve()