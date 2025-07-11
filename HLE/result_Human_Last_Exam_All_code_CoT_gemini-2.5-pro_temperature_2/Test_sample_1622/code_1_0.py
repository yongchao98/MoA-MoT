import sympy

def generate_formula():
    """
    Generates and prints the mathematical formula for P(n).
    """
    L, n = sympy.symbols('L n')

    # Coefficients for the n^-2 term derived from matching the expansions
    num2_L2 = 3
    num2_L1 = -2
    num2_L0 = 2
    den2 = 24
    term2 = (num2_L2*L**2 + num2_L1*L + num2_L0) / (den2*n**2)

    # Coefficients for the n^-3 term
    num3_L3 = 1
    num3_L2 = -2
    num3_L1 = 2
    den3 = 48
    term3 = (num3_L3*L**3 + num3_L2*L**2 + num3_L1*L) / (den3*n**3)

    # The full formula for P(n)
    P_n_formula = term2 + term3

    # Using sympy.pretty_print for a clear mathematical representation
    print("The derived formula for P(n) is:")
    # Sympy's pretty print can be wide, but it's the most mathematically correct.
    # For a single line output, we format it ourselves.
    
    # We will build a single-line string for the final answer block.
    # Note: Sympy's `str` output is used for its reliable formatting.
    str_term2 = f"({str(term2.args[1])}) / ({str(term2.args[0].args[0])})"
    str_term3 = f"({str(term3.args[1])}) / ({str(term3.args[0].args[0])})"
    
    # Python code output
    print(f"P(n) = {sympy.pretty(P_n_formula, use_unicode=False)}")
    
    # Final answer for capture
    final_answer_str = f"P(n) = {str(P_n_formula)}"
    
    # Optional: manual formatting for the desired output style
    manual_str = f"({num2_L2}*L^2 - {abs(num2_L1)}*L + {num2_L0}) / ({den2}*n^2) + ({num3_L3}*L^3 - {abs(num3_L2)}*L^2 + {num3_L1}*L) / ({den3}*n^3)"
    # Clean up 1*L^3 to L^3
    manual_str = manual_str.replace("1*L^3", "L^3")
    
    return manual_str


final_formula = generate_formula()

# The final answer must be wrapped in <<<>>>
print(f"\n<<<{final_formula}>>>")
