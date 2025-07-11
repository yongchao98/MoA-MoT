def solve():
    """
    Calculates the number of internal adjunctions in the given 2-category
    from F_11^3 to itself. This number is equivalent to the order of the
    general linear group GL(3, F_11).
    """
    n = 3
    q = 11

    # The order of GL(n, q) is (q^n - 1)(q^n - q)...(q^n - q^(n-1))
    terms = []
    q_n = q**n
    for i in range(n):
        term = q_n - q**i
        terms.append(term)

    result = 1
    for term in terms:
        result *= term

    print("The number of internal adjunctions is the order of the general linear group GL(3, F_11).")
    print("This is calculated as the product of the following numbers:")
    
    equation_str = " * ".join(map(str, terms))
    print(f"{equation_str} = {result}")

solve()