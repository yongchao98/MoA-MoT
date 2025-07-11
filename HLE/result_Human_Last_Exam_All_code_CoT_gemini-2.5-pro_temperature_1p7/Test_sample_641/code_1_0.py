def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The number of involutions in PSU(4, q) is given by the formula:
    # N = (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2
    # We calculate each term of the equation.

    # Calculate q^4
    term_q4 = q**4

    # Calculate q^2 - q + 1
    term_q2_q_1 = q**2 - q + 1

    # Calculate q^2 + 1
    term_q2_1 = q**2 + 1

    # The numerator of the formula is the number of non-central
    # involutions in SU(4, q).
    numerator = term_q4 * term_q2_q_1 * term_q2_1

    # The number of involutions in PSU(4, q) is half of that.
    num_involutions = numerator // 2

    # Output the equation with the computed numbers.
    print("The number of involutions is calculated by the equation:")
    print(f"N = (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2")
    print(f"For q = {q}:")
    print(f"N = ({term_q4} * {term_q2_q_1} * {term_q2_1}) / 2")
    print(f"N = {num_involutions}")


solve()