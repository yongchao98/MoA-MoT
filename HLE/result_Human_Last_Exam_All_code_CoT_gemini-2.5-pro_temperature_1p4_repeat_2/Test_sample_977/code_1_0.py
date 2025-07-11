def solve_potential():
    """
    This function prints the derived formula for the electric potential
    in the region 0 <= y <= a.
    """
    
    # Define the symbolic components of the formula as strings
    numerator = "-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    denominator = "k * [epsilon_2 * cosh(k*a) * sinh(k*b) + epsilon_1 * sinh(k*a) * cosh(k*b)]"

    print("The electric potential Phi(x, y) in the region 0 <= y <= a is given by:")
    print("")
    print("Phi(x, y) = " + numerator)
    print("            " + "-" * len(denominator))
    print("            " + denominator)

solve_potential()
