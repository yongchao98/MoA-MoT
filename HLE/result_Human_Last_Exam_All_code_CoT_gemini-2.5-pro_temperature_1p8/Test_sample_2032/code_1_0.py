import sympy

def solve_variance():
    """
    This function calculates the variance of Y using symbolic integration.
    Y is the second closest value to X1 among {X2, X3, X4}, where all Xi are i.i.d U(0,1).
    """

    # Define symbolic variables for x (a value for X1) and d (distance).
    x, d = sympy.symbols('x d')

    # Step 1: Define the PDF of the second order statistic of distances, conditional on X1=x.
    # We analyze for x in [0, 1/2] due to symmetry.
    
    # Part 1: For distances d in [0, x].
    # The CDF of the distance D=|X-x| is F(d) = 2d, and PDF is f(d)=2.
    F_D1 = 2 * d
    f_D1 = 2
    # The PDF of the 2nd order statistic from 3 samples is 6*F*(1-F)*f.
    pdf_D2_1 = 6 * F_D1 * (1 - F_D1) * f_D1

    # Part 2: For distances d in [x, 1-x].
    # The CDF of D is F(d) = x+d, and PDF is f(d)=1.
    F_D2 = x + d
    f_D2 = 1
    pdf_D2_2 = 6 * F_D2 * (1 - F_D2) * f_D2

    # Step 2: Define the conditional expectation of Y^2 given X1=x and D(2)=d.
    # Part 1: For d in [0, x], Y can be x+d or x-d with equal probability.
    # E[Y^2 | d] = 0.5*(x+d)^2 + 0.5*(x-d)^2 = x^2 + d^2
    integrand1 = (x**2 + d**2) * pdf_D2_1

    # Part 2: For d in [x, 1-x], Y must be x+d (since x-d < 0).
    # E[Y^2 | d] = (x+d)^2
    integrand2 = (x + d)**2 * pdf_D2_2

    # Step 3: Compute E[Y^2 | X1=x] by integrating over d. Let's call it G(x).
    G_integral_part1 = sympy.integrate(integrand1, (d, 0, x))
    G_integral_part2 = sympy.integrate(integrand2, (d, x, 1 - x))
    G_x = G_integral_part1 + G_integral_part2
    
    # Step 4: Compute the final unconditional expectation E[Y^2].
    # Due to symmetry G(x) = G(1-x), we integrate from 0 to 1/2 and multiply by 2.
    E_Y2 = 2 * sympy.integrate(G_x, (x, 0, sympy.S.Half))

    # Step 5: Calculate the variance using Var(Y) = E[Y^2] - (E[Y])^2.
    # We found E[Y] = 1/2 from a symmetry argument.
    E_Y = sympy.S.Half
    variance = E_Y2 - E_Y**2

    # Print the results in the required format.
    print("Step 1: Determine E[Y]")
    print("By symmetry of the U(0,1) distribution, E[Y] = 1/2.")
    print("")
    print("Step 2: Determine E[Y^2]")
    print(f"The conditional expectation E[Y^2|X1=x] for x in [0, 1/2] is: {sympy.simplify(G_x)}")
    print(f"Integrating this from x=0 to 1/2 and multiplying by 2 gives E[Y^2] = {E_Y2}")
    print("")
    print("Step 3: Calculate the variance")
    print(f"Var(Y) = E[Y^2] - (E[Y])^2")
    print(f"       = {E_Y2} - ({E_Y})^2")
    print(f"       = {E_Y2} - {E_Y**2}")
    print(f"       = {variance}")

if __name__ == '__main__':
    solve_variance()