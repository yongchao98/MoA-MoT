def find_blowup_conditions():
    """
    This function states the condition on the initial value y(0) for a given x(0) > 1
    that leads to a blow-up solution for the given system of differential equations.

    The final condition is derived analytically and presented as an inequality for y(0).
    The numbers in the equation are c1=2, c2=1, c3=3, p1=2, p2=3.
    """

    # The coefficients and powers in the expression for the bound of y(0)^2
    # y(0)^2 < c1*x(0) + c2 - c3*x(0)^(p1/p2)
    c1 = 2
    c2 = 1
    c3 = 3
    p1 = 2
    p2 = 3

    # Construct the string for the expression inside the square root
    bound_expression = f"{c1}*x(0) + {c2} - {c3}*x(0)^({p1}/{p2})"

    # Construct the final inequality string for y(0)
    inequality_string = f"-sqrt({bound_expression}) < y(0) < sqrt({bound_expression})"

    print("For a given x(0) > 1, the solution of the system blows up if y(0) satisfies the following inequality:")
    print(inequality_string)

find_blowup_conditions()