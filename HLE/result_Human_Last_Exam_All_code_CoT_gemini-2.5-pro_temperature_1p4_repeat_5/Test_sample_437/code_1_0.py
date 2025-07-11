import sympy

def solve_integral():
    """
    Symbolically evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.
    """
    # Define the symbolic variables for radial distance 'r' and orbital exponent 'zeta'.
    r, zeta = sympy.symbols('r zeta', real=True, positive=True)

    # The integral can be broken down into four parts based on the simplified expression:
    # (zeta**3/pi) * (Integral d(phi)) * (Integral sin(theta)d(theta)) * (Integral r*exp(-2*zeta*r)dr)

    # Part 1: Constant from the squared normalization of the 1s Slater orbital.
    norm_constant_part = sympy.sympify("zeta**3 / pi")

    # Part 2: The integral over the azimuthal angle 'phi' from 0 to 2*pi.
    phi_integral_val = sympy.integrate(1, (sympy.Symbol('phi'), 0, 2 * sympy.pi))

    # Part 3: The integral over the polar angle 'theta' from 0 to pi.
    theta_integral_val = sympy.integrate(sympy.sin(sympy.Symbol('theta')), (sympy.Symbol('theta'), 0, sympy.pi))

    # Part 4: The radial integral.
    radial_integrand = r * sympy.exp(-2 * zeta * r)
    radial_integral_val = sympy.integrate(radial_integrand, (r, 0, sympy.oo))

    # Now, we print out each component of the final equation.
    print("The integral is evaluated as a product of four terms:")
    print(f"1. Normalization constants squared: {norm_constant_part}")
    print(f"2. Integral over phi (azimuthal angle): {phi_integral_val}")
    print(f"3. Integral over theta (polar angle): {theta_integral_val}")
    print(f"4. Radial integral of (r * exp(-2*zeta*r)): {radial_integral_val}\n")

    # Calculate the final result by multiplying all the parts.
    final_result = norm_constant_part * phi_integral_val * theta_integral_val * radial_integral_val

    # Display the full equation with the value of each part.
    print("The final equation is the product of these terms:")
    final_equation = f"({norm_constant_part}) * ({phi_integral_val}) * ({theta_integral_val}) * ({radial_integral_val}) = {final_result}"
    print(final_equation)

    print("\nTherefore, the evaluated integral <phi_1s| 1/r |phi_1s> is:")
    print(final_result)

    print("\nNote: For a hydrogen atom in its ground state, the orbital exponent zeta = 1.")
    print(f"In this specific case, the value of the integral would be: {final_result.subs(zeta, 1)}")

solve_integral()