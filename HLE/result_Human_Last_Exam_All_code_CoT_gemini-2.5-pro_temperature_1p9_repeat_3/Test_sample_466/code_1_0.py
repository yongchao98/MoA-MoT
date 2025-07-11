def solve():
    """
    This function calculates the number of internal adjunctions in the specified
    2-category from F_11^3 to itself.
    This number is equivalent to the order of the general linear group GL(3, 11).
    """
    n = 3
    q = 11

    # The formula for the order of GL(n, F_q) is:
    # (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^{n-1})
    q_to_n = q**n

    term1_val = q_to_n - q**0
    term2_val = q_to_n - q**1
    term3_val = q_to_n - q**2

    result = term1_val * term2_val * term3_val

    # Output the components of the calculation as requested
    print(f"The number of adjunctions is the size of GL(3, F_11).")
    print(f"|GL(3, 11)| = ({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2)")
    print(f"             = ({q_to_n} - {q**0}) * ({q_to_n} - {q**1}) * ({q_to_n} - {q**2})")
    print(f"             = {term1_val} * {term2_val} * {term3_val}")
    print(f"             = {result}")

solve()