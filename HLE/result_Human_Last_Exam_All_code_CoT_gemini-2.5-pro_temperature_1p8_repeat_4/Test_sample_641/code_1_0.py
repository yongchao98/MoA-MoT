def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    The formula is (q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2.
    """
    q = 997

    # Calculate the components of the formula
    q_pow_4 = q**4
    q_sq_plus_1 = q**2 + 1
    q_sq_minus_q_plus_1 = q**2 - q + 1

    # Calculate the total number of involutions
    numerator = q_pow_4 * q_sq_plus_1 * q_sq_minus_q_plus_1
    result = numerator // 2

    # Print the equation with the calculated numbers and the final result
    print(f"Based on the formula N = (q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2:")
    print(f"For q = 997:")
    print(f"The equation is: ({q_pow_4} * {q_sq_plus_1} * {q_sq_minus_q_plus_1}) / 2")
    print(f"Number of involutions = {result}")

solve()