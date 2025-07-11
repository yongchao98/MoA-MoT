import sympy as sy

def solve_variance():
    """
    This function calculates the variance of Y symbolically using sympy.
    """
    # Define symbolic variables for y and x
    y, x = sy.symbols('y x', real=True, positive=True)

    # The PDF of Y, f_Y(y), is given by 6 * g(y), where g(y) is an integral.
    # We calculate g(y) by splitting the integral into 4 parts based on the piecewise
    # definition of the distance CDF.

    # Part 1: x in [0, y/2]
    integrand1 = y * (1 - y)
    g1 = sy.integrate(integrand1, (x, 0, y / 2))

    # Part 2: x in [y/2, y]
    integrand2 = (2 * y - 2 * x) * (1 - (2 * y - 2 * x))
    g2 = sy.integrate(integrand2, (x, y / 2, y))

    # Part 3: x in [y, (1+y)/2]
    integrand3 = (2 * x - 2 * y) * (1 - (2 * x - 2 * y))
    g3 = sy.integrate(integrand3, (x, y, (1 + y) / 2))

    # Part 4: x in [(1+y)/2, 1]
    integrand4 = y * (1 - y)
    g4 = sy.integrate(integrand4, (x, (1 + y) / 2, 1))

    # Combine the parts to get g(y)
    g_y = sy.simplify(g1 + g2 + g3 + g4)

    # The PDF of Y is f_Y(y) = 6 * g(y)
    f_y = sy.simplify(6 * g_y)

    # Calculate E[Y^2] by integrating y^2 * f_Y(y) from 0 to 1
    E_Y_sq = sy.integrate(y**2 * f_y, (y, 0, 1))

    # E[Y] is known to be 1/2 from symmetry
    E_Y = sy.S(1) / 2

    # Calculate the variance: Var(Y) = E[Y^2] - (E[Y])^2
    Var_Y = E_Y_sq - E_Y**2
    
    # Print the steps of the final calculation
    print("The variance of Y is calculated as Var(Y) = E[Y^2] - (E[Y])^2.")
    print(f"E[Y] = {E_Y}")
    print(f"E[Y^2] = {E_Y_sq}")
    print(f"Var(Y) = {E_Y_sq} - ({E_Y})^2 = {E_Y_sq} - {E_Y**2} = {Var_Y}")

solve_variance()