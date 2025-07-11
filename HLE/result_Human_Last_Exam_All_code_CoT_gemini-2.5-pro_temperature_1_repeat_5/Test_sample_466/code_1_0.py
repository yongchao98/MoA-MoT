def solve():
    """
    Calculates the number of internal adjunctions, which corresponds to the
    order of the general linear group GL(3, F_11).
    """
    q = 11
    n = 3

    # Calculate q^n
    q_n = q**n

    # The terms in the product formula for |GL(n, q)|
    term1 = q_n - (q**0)
    term2 = q_n - (q**1)
    term3 = q_n - (q**2)

    # Calculate the final result
    result = term1 * term2 * term3

    # Print the calculation steps
    print("An internal adjunction from V to itself is determined by an invertible linear map f: V -> V.")
    print("Thus, we need to count the number of such maps, which is the order of GL(3, F_11).")
    print("The formula for the order of GL(n, q) is: Product_{k=0 to n-1} (q^n - q^k)")
    print(f"For n=3 and q=11:")
    print(f"|GL(3, F_11)| = (11^3 - 11^0) * (11^3 - 11^1) * (11^3 - 11^2)")
    print(f"             = ({q_n} - {q**0}) * ({q_n} - {q**1}) * ({q_n} - {q**2})")
    print(f"             = {term1} * {term2} * {term3}")
    print(f"             = {result}")

solve()