def solve_order():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    n = 4
    q = 2

    # Step 1: Explain the theoretical background
    print("The defect group D is elementary abelian of order 16 (2^4).")
    print("This means D is isomorphic to a 4-dimensional vector space over the field F_2.")
    print("The inertial quotient E is a subgroup of the automorphism group Aut(D).")
    print("Therefore, E is a subgroup of the general linear group GL(4, F_2).")
    print("\nA key theorem in block theory states that since the field characteristic is 2, the order of E must be odd.")
    print("So, we need to find the largest odd divisor of the order of GL(4, 2).\n")

    # Step 2: Calculate the order of GL(4, 2)
    print("The formula for the order of GL(n, q) is (q^n - 1) * (q^n - q) * ... * (q^n - q^(n-1)).")
    
    terms = []
    for i in range(n):
        term = q**n - q**i
        terms.append(term)
    
    order_gl = 1
    for term in terms:
        order_gl *= term

    print(f"For GL(4, 2), the order is (2^4 - 2^0) * (2^4 - 2^1) * (2^4 - 2^2) * (2^4 - 2^3)")
    print(f"= ({terms[0]}) * ({terms[1]}) * ({terms[2]}) * ({terms[3]})")
    print(f"= {order_gl}\n")

    # Step 3: Find the largest odd divisor of the order
    print("To find the largest odd divisor, we can factorize the terms:")
    print(f"15 = 3 * 5")
    print(f"14 = 2 * 7")
    print(f"12 = 3 * 4 = 3 * 2^2")
    print(f"8 = 2^3")
    print("The odd factors are 3, 5, 3, and 7.")
    
    highest_order = order_gl
    while highest_order % 2 == 0:
        highest_order //= 2
    
    print("\nThe product of the odd factors gives the highest possible order for E.")
    # The final equation with each number as requested
    print(f"The highest possible order is 3^2 * 5 * 7 = 9 * 5 * 7 = {highest_order}")

solve_order()