def solve_work_equation():
    """
    This function programmatically constructs the formula for the work done
    by the current source based on the physical derivation.
    """

    # Symbolic components of the formula as strings
    term1 = "-(mu - mu_0)"
    term2 = "N**2"
    term3 = "w"
    term4 = "(x_2 - x_1)"
    term5 = "(I_2**2 - I_1**2)"
    denominator = "2*g"

    # Assemble the numerator and the full expression
    numerator_str = f"{term1} * {term2} * {term3} * {term4} * {term5}"
    work_formula_str = f"W = ({numerator_str}) / ({denominator})"

    # Let's print the final result in a format that closely matches option D.
    # W = - \frac{\mu - \mu_0}{2g} N^2 w (x_2 - x_1) (I_2^2 - I_1^2)
    # The code below will print out each symbolic component of the final equation.
    print("W = - ", end="")
    print("(mu - mu_0)", end="")
    print(" / ", end="")
    print("(2*g)", end="")
    print(" * ", end="")
    print("N**2", end="")
    print(" * ", end="")
    print("w", end="")
    print(" * ", end="")
    print("(x_2 - x_1)", end="")
    print(" * ", end="")
    print("(I_2**2 - I_1**2)", end="")
    print("")

solve_work_equation()