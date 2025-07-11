def solve_differential_equation():
    """
    This function prints the general solution of the given differential equation.
    The solution is derived by transforming the equation into a Clairaut equation
    using the substitution X=x**2 and Y=y**2.
    """
    
    # The coefficients and exponents in the final solution y**2 = C*x**2 + C**2 + 1*C + 9
    y_power = 2
    x_power = 2
    C_power = 2
    C_coefficient = 1
    constant_term = 9

    # Format the final equation string
    solution = f"y**{y_power} = C*x**{x_power} + C**{C_power} + {C_coefficient}*C + {constant_term}"

    print("The general solution of the differential equation is:")
    print(solution)

solve_differential_equation()
