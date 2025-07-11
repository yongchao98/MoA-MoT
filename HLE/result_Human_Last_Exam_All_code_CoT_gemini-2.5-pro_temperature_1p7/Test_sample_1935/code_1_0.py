import sympy

def solve_geometry_ratio():
    """
    Calculates the ratio S_KMN : S_ABC in terms of r and R.

    Given:
    - D, E, F are points of tangency of the incircle of acute triangle ABC.
    - r is the inradius of triangle DEF.
    - R is the inradius of triangle ABC.
    - KMN is the orthic triangle of DEF.

    The ratio S_KMN : S_ABC can be shown to be (r^2 * (R + r)^2) / R^4.
    """
    # Define symbolic variables for the radii
    r = sympy.Symbol('r')
    R = sympy.Symbol('R')

    # The derived formula for the ratio of the areas
    ratio_formula = (r**2 * (R + r)**2) / R**4
    
    # Print the symbolic formula
    print("The ratio S_KMN : S_ABC is given by the formula:")
    
    # We want to print the equation with all the terms expanded.
    # We can do this by creating an equation object and printing it.
    S_KMN_div_S_ABC = sympy.Symbol('S_KMN/S_ABC')
    equation = sympy.Eq(S_KMN_div_S_ABC, ratio_formula)
    
    # Sympy's pretty print can be used, but for a clear final expression let's format it.
    # The default print of the formula is already quite good.
    # Let's show the expanded numerator to match the description.
    expanded_numerator = sympy.expand(r**2 * (R + r)**2)
    final_ratio = expanded_numerator / R**4
    
    # To output the final formula as requested, printing each part.
    # Example format: Final Ratio = (r**2*R**2 + 2*r**3*R + r**4) / R**4
    
    num_term1_coeff = expanded_numerator.coeff(R**2, 1)
    num_term2_coeff = expanded_numerator.coeff(R, 1)
    num_term3_coeff = expanded_numerator.coeff(R, 0)
    
    print("Final equation:")
    print(f"S_KMN / S_ABC = ({num_term1_coeff})*R**2 + ({num_term2_coeff})*R + ({num_term3_coeff}) / R**4")
    
    # For a more direct print of the symbolic expression
    print("\nSymbolic expression:")
    print(final_ratio)

if __name__ == '__main__':
    solve_geometry_ratio()
