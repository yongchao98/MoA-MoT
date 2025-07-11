import sympy

def solve():
    """
    Computes the total Chern class of the Schur bundle S^{(2,1)}T(CP^4)
    as a polynomial in h.
    """
    # Setup the symbolic variable and the ring Q[h]/h^5
    h = sympy.Symbol('h')

    def truncate(poly):
        """Truncates a polynomial to degree 4 in h."""
        return sympy.series(poly, h, 0, 5).removeO()

    # Step 1: Find Chern classes of T = T(CP^4)
    # c(T) = (1+h)^5
    c_T_total = truncate((1 + h)**5)
    c_T_poly = sympy.Poly(c_T_total, h)
    
    # Extract Chern classes c_k(T)
    c = [c_T_poly.coeff_monomial(h**k) for k in range(5)]
    # c[0] is 1, c[k] corresponds to c_k for k>0
    c1, c2, c3, c4 = c[1], c[2], c[3], c[4]

    # Step 2: Compute Chern character ch(T)
    # The rank of T(CP^4) is 4
    rank_T = 4
    ch_T = [sympy.Integer(rank_T)]
    
    ch_T.append(c1)
    ch_T.append(sympy.Rational(1, 2) * (c1**2 - 2*c2))
    ch_T.append(sympy.Rational(1, 6) * (c1**3 - 3*c1*c2 + 3*c3))
    ch_T.append(sympy.Rational(1, 24) * (c1**4 - 6*c1**2*c2 + 3*c2**2 + 8*c1*c3 - 24*c4))
    
    ch_T_total = sum(ch_T[k] * h**k for k in range(5))

    # Step 3: Compute ch(S^{(2,1)}T) using the character formula
    # Let F = S^{(2,1)}T. ch(F) = 1/3 * (ch(T)^3 - ch(psi^3 T))
    
    # ch(psi^3 T) = sum(3^k * ch_k(T))
    ch_psi3_T = sum(3**k * ch_T[k] * h**k for k in range(5))
    
    # ch(T)^3
    ch_T_cubed = truncate(ch_T_total**3)
    
    # ch(F)
    ch_F_total = truncate(sympy.Rational(1, 3) * (ch_T_cubed - ch_psi3_T))
    
    ch_F_poly = sympy.Poly(ch_F_total, h)
    ch_F = [ch_F_poly.coeff_monomial(h**k) for k in range(5)]

    # Step 4: Convert ch(F) back to c(F)
    cf = [sympy.Integer(1)] # c0 is 1 for the total class
    
    # c_1(F)
    cf1 = ch_F[1]
    cf.append(cf1)
    
    # c_2(F)
    cf2 = sympy.Rational(1, 2) * (cf1**2 - 2*ch_F[2])
    cf.append(cf2)
    
    # c_3(F)
    cf3 = sympy.Rational(1, 3) * (6*ch_F[3] - cf1**3 + 3*cf1*cf2)
    cf.append(cf3)
    
    # c_4(F)
    P4 = cf1**4 - 6*cf1**2*cf2 + 3*cf2**2 + 8*cf1*cf3
    cf4 = sympy.Rational(1, 24) * (P4 - 24*ch_F[4])
    cf.append(cf4)

    # Step 5: Format and print the result
    c_F_total = cf[0] + cf[1]*h + cf[2]*h**2 + cf[3]*h**3 + cf[4]*h**4
    
    # Build the output string to match mathematical notation
    expr_str = str(cf[0])
    for i in range(1, 5):
        coeff = cf[i]
        if coeff == 0:
            continue
        
        sign = " + " if coeff > 0 else " - "
        coeff_abs = abs(coeff)
        
        if i == 1:
            h_part = "h"
        else:
            h_part = f"h^{i}"

        if coeff_abs == 1:
            term = f"{h_part}"
        else:
            term = f"{coeff_abs}*{h_part}"
        
        expr_str += sign + term

    print("The total Chern class c(S^{(2,1)}T(CP^4)) is expressed as a polynomial in h:")
    print(expr_str)

solve()