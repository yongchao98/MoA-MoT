def solve():
    """
    Solves the problem by calculating the value of a_{n,k,l} mod p based on the given recurrence relation
    and base-p representations of n, k, and l.
    """
    p = 21023

    # The recurrence relation is a_n,k,l = 12*a_{n-1,k,l} + 3*a_{n-1,k-1,l} + 75*a_{n-1,k,l-1} + 27*a_{n-1,k-2,l-2}.
    # The generating function for the sequence a_n is G_n(x,y) = P(x,y)^n, where
    # P(x,y) = 12 + 3*x + 75*y + 27*x^2*y^2.
    # a_{n,k,l} is the coefficient of x^k y^l in the expansion of P(x,y)^n.

    # n, k, and l are defined by their base-p digits (n_j, k_j, l_j), which are periodic with period 3:
    # For j = 0 mod 3: (n_j, k_j, l_j) = (5, 2, 2)
    # For j = 1 mod 3: (n_j, k_j, l_j) = (3, 1, 2)
    # For j = 2 mod 3: (n_j, k_j, l_j) = (2, 1, 1)

    # By a generalization of Lucas's Theorem, a_{n,k,l} mod p is the product of coefficients
    # C_j = [x^{k_j} y^{l_j}] P(x,y)^{n_j} mod p for each digit j.
    # The product consists of repeating factors C_0, C_1, C_2.

    # Step 1: Calculate the coefficients C_0, C_1, C_2 mod p.

    # For j=2 mod 3, we calculate C_2 = [x^1 y^1] P(x,y)^2 mod p.
    # The term x^1*y^1 in (12 + 3x + 75y + ...)^2 comes from 2 * (3x) * (75y).
    C_2 = (2 * 3 * 75) % p

    # For j=1 mod 3, we calculate C_1 = [x^1 y^2] P(x,y)^3 mod p.
    # The term x^1*y^2 comes from choosing 3x once and 75y twice. The multinomial coefficient is 3.
    # The coefficient is 3 * (3^1) * (75^2).
    C_1 = (3 * 3 * pow(75, 2, p)) % p

    # For j=0 mod 3, we calculate C_0 = [x^2 y^2] P(x,y)^5 mod p.
    # This coefficient arises from two combinations:
    # 1. Choosing 27x^2y^2 once and 12 four times. Multinomial coeff: 5!/(1!4!) = 5.
    #    Term: 5 * (27^1) * (12^4).
    # 2. Choosing 3x twice, 75y twice, and 12 once. Multinomial coeff: 5!/(2!2!1!) = 30.
    #    Term: 30 * (12^1) * (3^2) * (75^2).
    term1_C0 = (5 * 27 * pow(12, 4, p)) % p
    term2_C0 = (30 * 12 * pow(3, 2, p) * pow(75, 2, p)) % p
    C_0 = (term1_C0 + term2_C0) % p

    # Step 2: The number of factors of each type (C_0, C_1, C_2) is N = (3p+1)/2.
    # The final result is (C_0 * C_1 * C_2)^N mod p.
    X = (C_0 * C_1 * C_2) % p

    # Step 3: Calculate the exponent N.
    N = (3 * p + 1) // 2

    # Step 4: Calculate the final result using modular exponentiation.
    # To optimize, we can reduce the exponent modulo (p-1) by Fermat's Little Theorem.
    exponent = N % (p - 1)
    if exponent == 0 and X != 0:
        exponent = p - 1
    
    result = pow(X, exponent, p)

    # Print the explanation and the final equation.
    print(f"The prime for the calculation is p = {p}.")
    print("The problem requires calculating a_{n,k,l} mod p.")
    print("This is solved by finding coefficients of a polynomial raised to powers corresponding to base-p digits.")
    print("\nThe three repeating base coefficients are:")
    print(f"C_0 = [x^2 y^2] (12 + 3x + 75y + 27x^2y^2)^5 mod {p} = {C_0}")
    print(f"C_1 = [x^1 y^2] (12 + 3x + 75y + 27x^2y^2)^3 mod {p} = {C_1}")
    print(f"C_2 = [x^1 y^1] (12 + 3x + 75y + 27x^2y^2)^2 mod {p} = {C_2}")
    
    print("\nThe product of these coefficients is:")
    print(f"X = ({C_0} * {C_1} * {C_2}) mod {p} = {X}")
    
    print("\nEach coefficient appears N times in the full product, where:")
    print(f"N = (3*p + 1)/2 = {N}")
    
    print("\nThe final value is a_{n,k,l} = X^N mod p. The calculation is:")
    print(f"{X}^{N} mod {p} = {result}")

solve()
<<<20112>>>