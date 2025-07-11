import sympy as sp

def solve_and_print_ratio():
    """
    This function uses symbolic mathematics to derive the expression for the ratio BM/MI
    in terms of triangle side lengths a, b, and c, and then prints the final result.
    """
    # Define symbols for side lengths and semi-perimeter. They are positive quantities.
    a, b, c, s = sp.symbols('a b c s', positive=True)

    # From the geometric derivation, the ratio BI/MI can be expressed using angles as:
    # BI/MI = 2 * sin(A/2) * sin(C/2) / sin(B/2).
    # We use the half-angle formulas to convert this to side lengths.
    sin_A_2_sq = (s - b) * (s - c) / (b * c)
    sin_B_2_sq = (s - a) * (s - c) / (a * c)
    sin_C_2_sq = (s - a) * (s - b) / (a * b)

    # The square of the ratio BI/MI is therefore:
    ratio_BI_MI_sq = 4 * sin_A_2_sq * sin_C_2_sq / sin_B_2_sq
    
    # Sympy can simplify this algebraic expression.
    ratio_BI_MI_sq_simplified = sp.simplify(ratio_BI_MI_sq)
    
    # Take the square root to get BI/MI. Since side lengths are positive, the ratio is positive.
    ratio_BI_MI = sp.sqrt(ratio_BI_MI_sq_simplified)
    
    # The final desired ratio is BM/MI = BI/MI + 1.
    ratio_BM_MI = ratio_BI_MI + 1
    
    # Substitute the semi-perimeter s = (a+b+c)/2 into the expression.
    final_ratio_expr = ratio_BM_MI.subs(s, (a + b + c) / 2)
    
    # Simplify the final expression to get the result in terms of a, b, c.
    final_ratio = sp.simplify(final_ratio_expr)

    # The problem asks to output each "number" (which we interpret as each variable)
    # in the final equation.
    
    # Deconstruct the final symbolic expression to identify its components.
    numer, denom = final_ratio.as_numer_denom()
    numer_terms = numer.as_ordered_terms()
    
    print("The final expression for the ratio BM / MI is:")
    # Print a user-friendly version of the formula.
    print(f"({numer_terms[0]} + {numer_terms[1]}) / {denom}")

    print("\nThe components of the final equation are:")
    print(f"Component 1 in numerator: {numer_terms[0]}")
    print(f"Component 2 in numerator: {numer_terms[1]}")
    print(f"Component in denominator: {denom}")
    print("\nwhere 'a', 'b', and 'c' are the lengths of the sides opposite to vertices A, B, and C, respectively.")

solve_and_print_ratio()