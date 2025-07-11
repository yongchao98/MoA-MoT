def generate_milp_constraints():
    """
    This function generates and prints the two additional inequalities
    required to make the MILP system an exact encoding of f(x).

    The variables are x, y, a, b. M is a large positive constant.
    The inequalities are derived by modeling the piecewise nature of the function f(x)
    for the cases where x < 1.
    """
    # Inequality 1: Enforces y >= x when a=0 and b=0.
    # The general form is y >= x - M*b - M*a
    inequality1 = "y >= x - M*a - M*b"

    # Inequality 2: Enforces y >= 0 when a=0 and b=1.
    # The general form is y >= -M*(1-b) - M*a which expands to y >= M*b - M*a - M
    inequality2 = "y >= M*b - M*a - M"

    # Print the final inequalities in the requested format "y ~ A, y ~ B"
    print(f"{inequality1}, {inequality2}")

generate_milp_constraints()