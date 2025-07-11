def solve_ode():
    """
    This function prints the general solution to the given differential equation.
    The solution is y^2 = C*x^2 + C^2 + C + 9, where C is an arbitrary constant.
    """

    # The numbers present in the final equation are the powers and the constant term.
    power_of_y = 2
    power_of_x = 2
    power_of_C = 2
    constant_term = 9

    # The general solution is expressed as an equation.
    # C represents an arbitrary constant.
    solution_string = f"y**{power_of_y} = C*x**{power_of_x} + C**{power_of_C} + C + {constant_term}"

    print("The general solution of the differential equation is:")
    print(solution_string)

solve_ode()