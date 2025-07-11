import sympy as sp

def solve_for_K():
    """
    This function calculates the transformed stress-energy tensor component T_thetatheta
    and solves for K based on the equation given in the problem.
    """
    # 1. Define symbolic variables
    # T_cal represents the script T symbol in the problem
    r, theta, a, omega, T_cal = sp.symbols('r theta a omega T_cal', real=True, positive=True)

    # 2. Define coordinate transformation from polar to Cartesian
    x = r * sp.cos(theta)
    y = r * sp.sin(theta)

    # 3. Define the given Cartesian tensor components (spatial part)
    # The common factor is T_cal * (a*omega)**2
    C = T_cal * (a * omega)**2
    T_xx = C * sp.sin(theta)**2
    T_xy = -C * sp.sin(theta) * sp.cos(theta)
    T_yy = C * sp.cos(theta)**2

    # 4. Calculate the partial derivatives for the transformation law
    dx_dtheta = sp.diff(x, theta)
    dy_dtheta = sp.diff(y, theta)

    # 5. Apply the tensor transformation law for T'_thetatheta
    # T'_thetatheta = (dx/dtheta)^2 * T_xx + 2*(dx/dtheta)*(dy/dtheta)*T_xy + (dy/dtheta)^2 * T_yy
    T_prime_thetatheta = (dx_dtheta**2 * T_xx +
                         2 * dx_dtheta * dy_dtheta * T_xy +
                         dy_dtheta**2 * T_yy)

    # 6. The problem is for a ring of radius 'a', so we substitute r=a
    T_prime_thetatheta_at_a = T_prime_thetatheta.subs(r, a)

    # 7. Simplify the expression for the transformed component
    simplified_T_thetatheta = sp.simplify(T_prime_thetatheta_at_a)

    # 8. Solve for K using the given relation: T'_thetatheta = a^2*sin^2(theta)*T_cal + K
    K_symbol = sp.symbols('K')
    given_form = a**2 * sp.sin(theta)**2 * T_cal + K_symbol
    
    # The equation is: simplified_T_thetatheta = given_form
    # So, K = simplified_T_thetatheta - a**2*sin^2(theta)*T_cal
    K_expr = sp.solve(sp.Eq(simplified_T_thetatheta, given_form), K_symbol)[0]
    
    # Print the steps and the final result
    print("The transformed tensor component T_thetatheta is:")
    print(f"T'_θθ = {simplified_T_thetatheta}")
    print("\nThe problem states the equation:")
    # Using python symbols for the printout
    lhs, rhs = simplified_T_thetatheta, given_form
    print(f"{lhs} = {rhs}")
    print("\nSolving this equation for K, we get:")
    print(f"K = {K_expr}")

    # As requested, outputting the components/numbers of the final equation
    # LHS = simplified_T_thetatheta
    # RHS_part1 = a**2 * sp.sin(theta)**2 * T_cal
    # RHS_part2 = K
    print("\nThe final equation is built from these terms:")
    print(f"Term 1 (calculated T'_θθ): {simplified_T_thetatheta}")
    print(f"Term 2 (from given form): {a**2 * sp.sin(theta)**2 * T_cal}")
    print(f"Equating them gives K = Term 1 - Term 2.")

solve_for_K()