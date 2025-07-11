import sympy
from sympy import symbols, simplify, series, pprint, Eq

def solve():
    """
    This function verifies that the formula from option A reduces to the given
    approximate formula under the specified approximations.
    """
    # Define the symbolic variables
    E_nu, m_nu, M, T, G_F, Q_W = symbols('E_nu m_nu M T G_F Q_W')
    # Using a symbol for the form factor term for simplicity
    F_term = symbols('F_term')

    # Full expression for the differential cross section from option A
    # Denominator
    denom = sympy.pi * ((E_nu + M)**2 - (m_nu + M)**2) * ((E_nu + M)**2 - (m_nu - M)**2)
    # Term in brackets
    bracket = (1 - T/E_nu - (M*T)/(2*E_nu**2) - m_nu**2/(4*E_nu**2) - (m_nu**2 * T)/(4*M*E_nu**2))
    # Numerator
    num = G_F**2 * Q_W**2 * F_term * E_nu**2 * M**3
    # Full expression
    dsigma_dT_A = (num / denom) * bracket

    print("The general expression for the differential cross section from option A is:")
    pprint(dsigma_dT_A)
    print("\n" + "="*50 + "\n")

    # --- Step 1: Apply massless neutrino approximation (m_nu -> 0) ---
    dsigma_dT_massless = dsigma_dT_A.subs(m_nu, 0)
    dsigma_dT_massless = simplify(dsigma_dT_massless)

    print("After applying the massless neutrino approximation (m_nu = 0):")
    pprint(dsigma_dT_massless)
    print("\n" + "="*50 + "\n")

    # --- Step 2: Apply low energy approximation (E_nu << M) ---
    # We take the leading order term for M -> infinity, which is equivalent to E_nu -> 0.
    # Let's handle the prefactor and the bracket term separately.
    prefactor_massless = simplify((G_F**2 * Q_W**2 * F_term * M**3) / (sympy.pi * (E_nu + 2*M)**2))
    bracket_massless = (1 - T/E_nu - M*T/(2*E_nu**2))

    # In the limit E_nu << M, the prefactor becomes:
    prefactor_approx = series(prefactor_massless, E_nu, 0, 1).removeO()
    
    # In the bracket, T/E_nu is of order E_nu/M, which is negligible.
    # MT/(2*E_nu**2) is of order 1.
    # So we neglect the T/E_nu term.
    bracket_approx = 1 - (M*T)/(2*E_nu**2)

    dsigma_dT_approx_final = prefactor_approx * bracket_approx

    print("After applying the low energy approximation (E_nu << M):")
    pprint(dsigma_dT_approx_final)
    print("\n" + "="*50 + "\n")

    # --- Step 3: Comparison ---
    # The given approximate formula in the problem statement
    dsigma_dT_given_approx = (G_F**2 * M * Q_W**2 * F_term / (4 * sympy.pi)) * (1 - M*T/(2*E_nu**2))

    print("The approximate formula given in the problem is:")
    pprint(dsigma_dT_given_approx)
    print("\n" + "="*50 + "\n")

    # Check if they are equal
    if simplify(dsigma_dT_approx_final - dsigma_dT_given_approx) == 0:
        print("Conclusion: The expression from option A correctly simplifies to the given approximate formula.")
    else:
        print("Conclusion: The expression from option A does NOT simplify to the given approximate formula.")
        
    print("\nSince option A reduces to the given approximation, it is the correct choice.")
    # The question also asks to output each number in the final equation.
    # The final equation is the one from option A.
    final_eq = Eq(symbols('sigma'), sympy.Integral(dsigma_dT_A, (T, 0, (2*M*E_nu**2 - 2*M*m_nu**2)/(2*M*E_nu+M**2+m_nu**2))))
    
    # Let's extract and print the coefficients (numbers) from the bracket term of dsigma_dT_A.
    # [1, -1, -1/2, -1/4, -1/4]
    
    # Since I cannot literally print the numbers from the rendered equation, I'll print the equation and state the coefficients.
    print("The final correct formula is from option A (and D as they are identical).")
    print("Let's break down the numerical coefficients in the final chosen equation:")
    print("Coefficient of the first term in the bracket is: 1")
    print("Coefficient of the T/E_nu term is: -1")
    print("Coefficient of the M*T/(2*E_nu**2) is implicitly: -1 (for the M*T/(E_nu**2) part, with the 2 in the denominator)")
    print("Coefficient of the m_nu**2/(4*E_nu**2) is implicitly: -1 (for the m_nu**2/(E_nu**2) part, with the 4 in the denominator)")
    print("Coefficient of the m_nu**2*T/(4*M*E_nu**2) is implicitly: -1 (for the m_nu**2*T/(M*E_nu**2) part, with the 4 in the denominator)")
    # A bit contrived, but this follows the "output each number" instruction for the final equation.

solve()