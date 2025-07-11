import sympy

def solve_and_print():
    """
    This function determines the coefficients c1 and c2 based on the reasoning
    that the bound should be tight for the simple case K=I.
    """
    c1 = 1
    c2 = 1

    # Using sympy to format the LaTeX output for the equation
    beta, alpha_i, K_alpha_i, K_u_i = sympy.symbols(
        r'\beta \alpha^{\mathcal{D}}_i (K\vec{\alpha}^{\mathcal{D}})_i -(K\vec{\alpha}^{\mathcal{D}-i})_i'
    )
    
    # Construct the inequality string
    lhs = sympy.latex(K_u_i)
    
    # Construct the RHS terms with the determined coefficients
    term1 = f"(1 + {c1} {sympy.latex(beta)}){sympy.latex(alpha_i)}"
    term2 = f"(1 + {c2} {sympy.latex(beta)}){sympy.latex(K_alpha_i)}"
    
    print(f"Based on the analysis, the coefficients are:")
    print(f"c_1 = {c1}")
    print(f"c_2 = {c2}")
    print("\nThe extended Jaakola-Haussler bound is:")
    
    # Print the full inequality using the derived coefficients
    # Python's f-strings require doubling the curly braces to escape them.
    # The LaTeX code also needs escaping for backslashes.
    
    final_equation = (
        f"-(K \\vec{{\\alpha}}^{{\\mathcal{{D}}-i}})_i "
        f"\\leq "
        f"+(1 + {c1}\\beta)\\alpha^{{\\mathcal{{D}}}}_i "
        f"- (1 + {c2}\\beta)(K \\vec{{\\alpha}}^{{\\mathcal{{D}}}})_i "
        f"+ o(\\beta)"
    )
    
    print(final_equation)

solve_and_print()