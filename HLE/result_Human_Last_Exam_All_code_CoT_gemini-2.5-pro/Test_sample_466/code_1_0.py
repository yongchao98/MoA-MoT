def solve():
    """
    Calculates the number of internal adjunctions in the specified 2-category
    from F_11^3 to itself.

    This number is equal to the order of the general linear group GL(3, 11).
    """
    q = 11
    n = 3

    # The order of GL(n, q) is (q^n - q^0)(q^n - q^1)...(q^n - q^(n-1))
    term1 = q**n - q**0
    term2 = q**n - q**1
    term3 = q**n - q**2

    result = term1 * term2 * term3

    print("The number of internal adjunctions is the order of the group GL(n, q), for n=3 and q=11.")
    print("The formula is |GL(3, 11)| = (11^3 - 11^0) * (11^3 - 11^1) * (11^3 - 11^2).")
    print(f"The first term in the equation is {q**n} - {q**0} = {term1}")
    print(f"The second term in the equation is {q**n} - {q**1} = {term2}")
    print(f"The third term in the equation is {q**n} - {q**2} = {term3}")
    print("The final equation is:")
    print(f"{term1} * {term2} * {term3} = {result}")

solve()