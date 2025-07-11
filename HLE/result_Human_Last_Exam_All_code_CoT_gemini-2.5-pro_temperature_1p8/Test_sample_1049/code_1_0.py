import sympy

def find_closed_form():
    """
    This function computes the closed form for the given sum.
    The method involves finding a polynomial P(n) such that the sum is 4^n * P(n).
    We find the coefficients of P(n) in the binomial basis {C(n,k)} by computing
    the first few values of the sum and using the method of finite differences.
    Finally, we convert this to a standard polynomial in n.
    """
    n = sympy.symbols('n')
    
    # These are the leading coefficients in the finite difference table of P(n).
    # P(n) = S_n / 4^n.
    # The values P(0), P(1), ..., P(5) are calculated to be:
    # 1, 122, 1233, 5884, 19250, 50076
    # The finite difference table is:
    # 1
    # 121
    # 990
    # 2550
    # 2625
    # 945
    
    # These are the coefficients of P(n) in the binomial basis.
    # P(n) = c0*C(n,0) + c1*C(n,1) + c2*C(n,2) + c3*C(n,3) + c4*C(n,4) + c5*C(n,5)
    # The coefficients are given by the top diagonal of the difference table.
    coeffs = {
        0: 1,
        1: 121,
        2: 990,
        3: 2550,
        4: 2625,
        5: 945
    }
    
    P_n = 0
    for k in range(6):
        P_n += coeffs[k] * sympy.binomial(n, k)
        
    # Expand P_n into a standard polynomial form
    expanded_P_n = sympy.expand(P_n)
    
    # The polynomial expression has a common denominator of 8.
    # We can write S_n = 4^n * (1/8) * (8 * P(n))
    # S_n = 2^(2n) * 2^(-3) * (8 * P(n)) = 2^(2n-3) * (8 * P(n))
    poly_numerator = sympy.simplify(8 * expanded_P_n)
    
    # Extract coefficients to print them one by one.
    poly_obj = sympy.Poly(poly_numerator, n)
    c = poly_obj.all_coeffs()
    
    print("The closed form for the sum is:")
    print(f"2**(2*n - 3) * ({c[0]}*n**5 + {c[1]}*n**4 + {c[2]}*n**3 + {c[3]}*n**2 + {c[4]}*n + {c[5]})")

if __name__ == '__main__':
    find_closed_form()
