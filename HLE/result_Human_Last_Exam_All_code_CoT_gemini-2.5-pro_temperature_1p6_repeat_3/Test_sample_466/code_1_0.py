import math

def solve():
    """
    Calculates the number of internal adjunctions by computing the order of GL(3, F_11).
    """
    n = 3
    q = 11

    # The number of such adjunctions is the number of invertible 3x3 matrices over F_11,
    # which is the order of the general linear group GL(3, 11).
    # The formula for |GL(n, q)| is (q^n - 1)(q^n - q)...(q^n - q^(n-1)).

    q_n = q**n

    term1_val = q_n - q**0
    term2_val = q_n - q**1
    term3_val = q_n - q**2

    result = term1_val * term2_val * term3_val

    print("The number of internal adjunctions is the order of the general linear group GL(n, q), for n=3 and q=11.")
    print("The formula is |GL(3, 11)| = (11^3 - 11^0) * (11^3 - 11^1) * (11^3 - 11^2).")
    print(f"This evaluates to:")
    print(f"({q_n} - {q**0}) * ({q_n} - {q**1}) * ({q_n} - {q**2})")
    print(f"= {term1_val} * {term2_val} * {term3_val}")
    print(f"= {result}")

solve()