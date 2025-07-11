def solve_inertial_quotient_order():
    """
    Calculates the highest possible order for the inertial quotient E of a block B
    with an elementary abelian defect group D of order 16 over a field of characteristic 2.
    """
    n = 4  # Dimension of the vector space, from D being (C_2)^4
    q = 2  # Size of the field, from the characteristic being 2

    print("The highest possible order for the inertial quotient E is the order of GL(n, q), where n=4 and q=2.")
    print("The formula for the order of GL(n, q) is:")
    print("|GL(n, q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))")
    print("\nFor n=4 and q=2, the calculation is:")

    # Calculate each term in the product
    term1 = q**n - q**0
    term2 = q**n - q**1
    term3 = q**n - q**2
    term4 = q**n - q**3
    
    # Calculate the intermediate powers
    q_n = q**n
    q_0 = q**0
    q_1 = q**1
    q_2 = q**2
    q_3 = q**3

    print(f"|GL(4, 2)| = (2^4 - 2^0) * (2^4 - 2^1) * (2^4 - 2^2) * (2^4 - 2^3)")
    print(f"           = ({q_n} - {q_0}) * ({q_n} - {q_1}) * ({q_n} - {q_2}) * ({q_n} - {q_3})")
    
    # The final equation with each number explicitly shown
    print(f"           = {term1} * {term2} * {term3} * {term4}")

    # Calculate the final result
    result = term1 * term2 * term3 * term4
    print(f"           = {result}")

    print("\nThus, the highest order that E can have is 20160.")

solve_inertial_quotient_order()