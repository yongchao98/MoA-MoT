import sympy

def find_maclaurin_coefficient():
    """
    This function calculates the 4th Maclaurin series coefficient of the given function
    using the sympy library for symbolic mathematics.
    """
    # Define the symbolic variable x and constants
    x = sympy.symbols('x')
    E = sympy.E
    pi = sympy.pi
    S = sympy.S

    # The function f(x) is the sum of two terms, T1(x) + T2(x).
    # The total coefficient is the sum of the coefficients from each term.

    # --- Term 1: T1(x) = 9 * x**4 / (16 * e) ---
    # The coefficient of x^4 in T1(x) is trivial.
    coeff_t1 = S(9) / (16 * E)

    # --- Term 2: The complex fraction T2(x) ---
    # Define the numerator and denominator parts of T2(x)
    
    # Numerator parts
    N1 = x**4 - S(5)/6 * sympy.log(x**4 + 1)**2
    N2 = sympy.exp(sympy.tanh(x**3)/2) - 1
    N3 = sympy.cos(sympy.sin(pi * sympy.cosh(x**6))) - 1/E
    Numerator = 4 * N1 * N2 * N3

    # Denominator parts
    D1 = sympy.tan(x**6) - sympy.log(x**8 + 1)
    D2 = sympy.exp(sympy.cos(x**5)**2 + sympy.sinh(x**2)) - 1
    D3 = sympy.cosh(x**3) - sympy.sec(x**7)
    Denominator = D1 * D2 * D3

    # The series for T2(x) is a Laurent series starting with O(x^-5). To find the
    # x^4 coefficient, we need to compute the series expansions to a high enough order.
    # Required order = (target_power) - (series_order) + (denominator_order)
    # = 4 - (-5) + 12 = 21. Let's use order 25 for safety.
    order = 25

    # Compute series for Numerator and Denominator separately
    num_series_poly = Numerator.series(x, 0, order).removeO()
    den_series_poly = Denominator.series(x, 0, order).removeO()

    # Create an approximate rational function from the series polynomials
    T2_approx = num_series_poly / den_series_poly

    # Compute the Laurent series of the T2 approximation. We need terms up to x^4.
    t2_series = T2_approx.series(x, 0, 5)

    # Extract the coefficient of x^4 from the series of T2(x)
    coeff_t2 = t2_series.coeff(x**4)

    # --- Final Result ---
    # The total coefficient is the sum of the coefficients from T1 and T2.
    total_coeff = coeff_t1 + coeff_t2
    
    # Print the components of the final calculation as requested
    print(f"The 4th Maclaurin coefficient is the sum of coefficients from two terms.")
    print(f"Coefficient from the first term: {sympy.pretty(coeff_t1)}")
    print(f"Coefficient from the second term: {coeff_t2}")
    print(f"Total Coefficient = {sympy.pretty(coeff_t1)} + {coeff_t2} = {sympy.pretty(total_coeff)}")

if __name__ == '__main__':
    find_maclaurin_coefficient()
