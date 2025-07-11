import sympy as sp

def find_maclaurin_coefficient():
    """
    Finds the 4th Maclaurin series coefficient of the given function by analyzing
    its components.
    """
    x = sp.Symbol('x')

    # --- Introduction and Part A ---
    print("### Step 1: Decompose the function")
    print("The function f(x) is split into two parts: f(x) = A(x) + B(x).")
    
    A_x = (9 * x**4) / (16 * sp.exp(1))
    print(f"\nA(x) = {sp.pretty(A_x, use_unicode=False)}")
    print("B(x) is the complex fractional part.")
    print("\nThe 4th Maclaurin coefficient of f(x) is the sum of the coefficients of x^4 from A(x) and B(x).")

    coeff_A = A_x.coeff(x**4)
    print("\n### Step 2: Analyze A(x)")
    print("A(x) is a monomial, so its x^4 coefficient is the constant multiplying x^4.")
    print(f"Coefficient of x^4 in A(x) = {sp.pretty(coeff_A, use_unicode=False)}")

    # --- Part B Analysis ---
    print("\n### Step 3: Analyze B(x) near x = 0")
    print("We determine the lowest power of x in the series expansion of B(x).")
    
    # Define components of B(x)
    # Numerator parts
    N1 = x**4 - sp.Rational(5, 6) * sp.log(x**4 + 1)**2
    N2 = sp.exp(sp.tanh(x**3) / 2) - 1
    N3 = sp.cos(sp.sin(sp.pi * sp.cosh(x**6))) - 1 / sp.exp(1)
    
    # Denominator parts
    D1 = sp.tan(x**6) - sp.log(x**8 + 1)
    D2 = sp.exp(sp.cos(x**5)**2 + sp.sinh(x**2)) - 1
    D3 = sp.cosh(x**3) - sp.sec(x**7)

    # Find the order of the leading term for each component
    def get_order(expr):
        leading_term = expr.as_leading_term(x)
        return sp.degree(leading_term, gen=x)

    order_N1 = get_order(N1)
    order_N2 = get_order(N2)
    order_N3 = get_order(N3)
    order_D1 = get_order(D1)
    order_D2 = get_order(D2)
    order_D3 = get_order(D3)
    
    print("\nAnalysis of the Numerator of B(x):")
    print(f"  - Lowest power in (x**4 - 5/6*log(x**4+1)**2) is x^{order_N1}.")
    print(f"  - Lowest power in (exp(tanh(x**3)/2) - 1) is x^{order_N2}.")
    print(f"  - Lowest power in (cos(sin(pi*cosh(x**6))) - 1/e) is x^{order_N3} (i.e., a non-zero constant).")
    order_num = order_N1 + order_N2 + order_N3
    print(f"  => Total lowest power in Numerator = {order_N1} + {order_N2} + {order_N3} = {order_num}.")

    print("\nAnalysis of the Denominator of B(x):")
    print(f"  - Lowest power in (tan(x**6) - log(x**8+1)) is x^{order_D1}.")
    print(f"  - Lowest power in (exp(cos(x**5)**2 + sinh(x**2)) - 1) is x^{order_D2} (i.e., a non-zero constant).")
    print(f"  - Lowest power in (cosh(x**3) - sec(x**7)) is x^{order_D3}.")
    order_den = order_D1 + order_D2 + order_D3
    print(f"  => Total lowest power in Denominator = {order_D1} + {order_D2} + {order_D3} = {order_den}.")
    
    order_B = order_num - order_den
    print(f"\nThe series for B(x) begins with a term of order x^({order_num} - {order_den}) = x^({order_B}).")
    print("Since the series for B(x) starts with a negative power of x, it is a Laurent series with a pole at x=0.")
    print("Such a series has no terms with non-negative powers, including x^4.")
    
    coeff_B = 0
    print(f"Coefficient of x^4 in B(x) = {coeff_B}")

    # --- Final Result ---
    print("\n### Step 4: Final Calculation")
    final_coeff = coeff_A + coeff_B
    
    pretty_A = sp.pretty(coeff_A, use_unicode=False)
    pretty_final = sp.pretty(final_coeff, use_unicode=False)

    print("The final coefficient is the sum of the coefficients from A(x) and B(x).")
    # Using 'pretty' formatting might create multiline output, so we format it manually
    num_A, den_A = sp.fraction(coeff_A)
    num_final, den_final = sp.fraction(final_coeff)
    
    print(f"\nFinal Equation:")
    print(f"4th Maclaurin coefficient = ({num_A}/{den_A}) + {coeff_B} = {num_final}/{den_final}")

if __name__ == '__main__':
    find_maclaurin_coefficient()
