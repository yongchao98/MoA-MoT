def solve():
    """
    Calculates the number of internal adjunctions from F_11^3 to itself
    in the given 2-category C.
    This corresponds to calculating the order of the general linear group GL(3, F_11).
    """
    n = 3
    q = 11

    # The order of GL(n, q) is (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))
    
    # Calculate each term in the product
    term1 = q**n - q**0
    term2 = q**n - q**1
    term3 = q**n - q**2
    
    # Calculate the final result
    result = term1 * term2 * term3
    
    print("The number of internal adjunctions is the order of the group GL(n, F_q).")
    print(f"The calculation for n={n} and q={q} is:")
    print(f"({q}^{n} - {q}^{0}) * ({q}^{n} - {q}^{1}) * ({q}^{n} - {q}^{2})")
    print(f"= ({q**n} - {q**0}) * ({q**n} - {q**1}) * ({q**n} - {q**2})")
    # Print the final equation with all numbers
    print(f"= {term1} * {term2} * {term3}")
    print(f"= {result}")

solve()