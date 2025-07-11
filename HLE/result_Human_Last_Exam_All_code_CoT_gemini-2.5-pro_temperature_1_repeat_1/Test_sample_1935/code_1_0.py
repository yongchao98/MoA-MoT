def solve_geometry_ratio():
    """
    This function prints the symbolic formula for the ratio S_KMN : S_ABC.
    """
    # The ratio of the area of triangle KMN to the area of triangle ABC
    # is found to be (r/R)^4, where r is the inradius of triangle DEF
    # and R is the inradius of triangle ABC.

    # Let's represent the formula as a string.
    numerator_base = "r"
    numerator_power = "4"
    denominator_base = "R"
    denominator_power = "4"

    print("The ratio S_KMN : S_ABC is given by the formula:")
    print(f"({numerator_base}/{denominator_base})^{numerator_power}")
    # An alternative way to print the expanded formula
    # print(f"{numerator_base}^{numerator_power} / {denominator_base}^{denominator_power}")

solve_geometry_ratio()