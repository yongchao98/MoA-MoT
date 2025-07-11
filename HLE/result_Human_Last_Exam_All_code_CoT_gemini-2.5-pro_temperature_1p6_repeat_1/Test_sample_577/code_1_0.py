import sympy

def compute_chern_class_S21_TCP4():
    """
    Computes the total Chern class of the bundle S^(2,1)T(CP^4).

    The plan is as follows:
    1. Define the generator h of the cohomology ring of CP^4.
    2. Compute the Chern character of the tangent bundle T(CP^4), denoted ch_V.
       This is derived from the Euler sequence: ch(T(CP^4)) = 5*exp(h) - 1.
    3. Compute the Chern characters of Adams operations psi^k(V), denoted ch_psik_V.
       This is done by scaling the j-th degree term of ch(V) by k^j.
    4. Use these to find the Chern characters of the exterior powers Lambda^2(V) and Lambda^3(V).
       ch(L^2 V) = 1/2 * (ch(V)^2 - ch(psi^2 V))
       ch(L^3 V) = 1/6 * (ch(V)^3 - 3*ch(V)*ch(psi^2 V) + 2*ch(psi^3 V))
    5. Find the Chern character of S^(2,1)V using the K-theory identity:
       S^(2,1)V = V tensor L^2(V) - L^3(V).
       ch(S^(2,1)V) = ch(V) * ch(L^2 V) - ch(L^3 V).
    6. Convert the resulting Chern character ch(E) to Chern classes c_k(E).
       The relations are given by Newton's sums.
       c1 = ch1
       c2 = (c1^2 - 2*ch2)/2
       c3 = (c1^3 - 3*c1*c2 + 6*ch3)/3
       c4 = (c1^4 - 6*c1^2*c2 + 8*c1*c3 + 3*c2^2 - 24*ch4)/6
    7. Print the total Chern class as 1 + c1 + c2 + c3 + c4.
    """
    h = sympy.Symbol('h')
    
    # The calculations are in the ring Z[h]/(h^5)
    def series(p, n=5):
        return sympy.series(p, h, 0, n).removeO()

    # Step 2: Chern character of V = T(CP^4)
    ch_V = series(5 * sympy.exp(h) - 1)
    ch_V_coeffs = [ch_V.coeff(h, i) for i in range(5)]
    
    # Step 3: Adams operations
    def ch_psi(p):
        res = 0
        for i, coeff in enumerate(ch_V_coeffs):
            res += coeff * (p**i) * h**i
        return res
        
    ch_psi2_V = ch_psi(2)
    ch_psi3_V = ch_psi(3)
    
    # Step 4: Exterior powers
    ch_L2_V = series(sympy.Rational(1, 2) * (ch_V**2 - ch_psi2_V))
    ch_L3_V = series(sympy.Rational(1, 6) * (ch_V**3 - 3*ch_V*ch_psi2_V + 2*ch_psi3_V))
    
    # Step 5: S^(2,1)V
    E = S21V = 'S^(2,1)T(CP^4)'
    ch_E = series(ch_V * ch_L2_V - ch_L3_V)
    ch = [ch_E.coeff(h, i) * h**i for i in range(5)]

    # Step 6: Convert ch(E) to c_k(E)
    c = [0] * 5
    c[1] = ch[1]
    # c2 from 2*ch2 = c1^2 - 2c2
    c[2] = series((c[1]**2 - 2*ch[2])/2)
    # c3 from 6*ch3 = c1^3 - 3*c1*c2 + 3*c3
    c[3] = series((c[1]**3 - 3*c[1]*c[2] + 6*ch[3])/3)
    # c4 from 24*ch4 = c1^4 - 6*c1^2*c2 + 8*c1*c3 + 3*c2^2 - 6*c4
    c[4] = series((c[1]**4 - 6*c[1]**2*c[2] + 8*c[1]*c[3] + 3*c[2]**2 - 24*ch[4])/6)

    # Step 7: Print the result
    total_c = 1 + c[1] + c[2] + c[3] + c[4]
    
    print(f"The total Chern class of {E} is:")
    print(f"c({E}) = 1 + ({sympy.poly(c[1], h).coeffs()[0]})h + ({sympy.poly(c[2], h).coeffs()[0]})h^2 + ({sympy.poly(c[3], h).coeffs()[0]})h^3 + ({sympy.poly(c[4], h).coeffs()[0]})h^4")

compute_chern_class_S21_TCP4()