def solve():
    """
    Calculates the number of internal adjunctions in the given 2-category from F_11^3 to itself.
    This number is the order of the general linear group GL(3, F_11).
    """
    q = 11
    n = 3

    q_n = q**n
    q_0 = q**0
    q_1 = q**1
    q_2 = q**2

    term1 = q_n - q_0
    term2 = q_n - q_1
    term3 = q_n - q_2

    result = term1 * term2 * term3

    print(f"({q}^{n} - {q}^{0}) * ({q}^{n} - {q}^{1}) * ({q}^{n} - {q}^{2}) = ({q_n} - {q_0}) * ({q_n} - {q_1}) * ({q_n} - {q_2}) = {term1} * {term2} * {term3} = {result}")

solve()