import sympy

def solve_for_chi_star():
    """
    Solves for the susceptibility chi* based on the duality principle in magnetostatics.

    The principle states that for a 2D system, a configuration with a material of
    susceptibility chi is dual to a configuration with susceptibility chi_star
    such that their relative permeabilities multiply to 1.
    
    (1 + chi) * (1 + chi_star) = 1
    
    This function solves this equation for chi_star.
    """
    # Define chi and chi_star as symbolic variables
    chi = sympy.Symbol('chi')
    chi_star = sympy.Symbol('chi_star')

    # Define the equation based on the duality principle
    # mu_r * mu_r_star = 1  => (1 + chi) * (1 + chi_star) = 1
    equation = sympy.Eq((1 + chi) * (1 + chi_star), 1)

    # Solve the equation for chi_star
    solution = sympy.solve(equation, chi_star)

    # The result is a list with one element, which is the expression for chi_star
    chi_star_expr = solution[0]
    
    # Print the final equation
    final_equation = sympy.Eq(chi_star, chi_star_expr)
    print("The relationship between chi* and chi is:")
    # Using sympy.pretty_print for a nicer output format
    sympy.pretty_print(final_equation)
    
    # The final expression is chi_star = -chi / (1 + chi)
    # The required format is <<<answer content>>>.
    # We output the final expression as a string.
    # The problem asks to output each number, so we construct the string carefully.
    numerator_str = f"-{chi}"
    denominator_str = f"{1} + {chi}"
    
    # To satisfy "output each number in the final equation", let's print the parts.
    print("\nFinal equation breakdown:")
    print(f"chi* = ({numerator_str}) / ({denominator_str})")


solve_for_chi_star()

# Final answer format for parsing
chi_symbol = sympy.Symbol('chi')
result_expr = -chi_symbol / (1 + chi_symbol)
print(f"<<<{result_expr}>>>")
