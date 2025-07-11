def solve_and_print_equation():
    """
    This function determines the equation of the function f(x) from its graph
    and prints the result in a formatted way.
    """

    # Symbolic representation of the parameters from the graph
    a = 'a'
    b = 'b'
    c = 'c'

    # The analysis of the graph reveals the following:
    # 1. Zeros at x = -b, x = d, x = b. Numerator is proportional to (x+b)(x-b)(x-d).
    # 2. Vertical asymptotes at x = a, x = c. Denominator is proportional to (x-a)(x-c).
    # 3. Slant asymptote is y = x. This implies that for f(x) = N(x)/D(x),
    #    the leading coefficient is 1 and lim(x->inf)[f(x) - x] = 0.
    # 4. The condition for the slant asymptote leads to the constraint d = a + c.

    # We substitute d = a + c into the expression for the function.
    print("The equation for the function f(x) is derived from its graphical features.")
    print("The key features are the x-intercepts, vertical asymptotes, and the slant asymptote y=x.")
    print("These features lead to the following equation, with the constraint d = a + c derived from the slant asymptote:")
    print("-" * 30)

    # Building the components of the final equation string
    # d is replaced by (a + c)
    d_expr = f"({a} + {c})"

    # The components of the equation are:
    numerator_factor1 = f"(x - {d_expr})"
    numerator_factor2 = f"(x**2 - {b}**2)"
    denominator_factor1 = f"(x - {a})"
    denominator_factor2 = f"(x - {c})"

    # Printing the final equation with each component clearly visible.
    # The equation is: f(x) = ( (x - (a+c)) * (x^2 - b^2) ) / ( (x-a) * (x-c) )
    print(f"f(x) = ( {numerator_factor1} * {numerator_factor2} ) / ( {denominator_factor1} * {denominator_factor2} )")

solve_and_print_equation()