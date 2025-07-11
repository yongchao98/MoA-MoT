import sympy

def evaluate_slater_integral():
    """
    Symbolically evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital
    using the sympy library.
    """

    # Define the mathematical symbols required for the calculation.
    # r, theta, phi are the spherical coordinates.
    # zeta is the orbital exponent, defined as a positive real number.
    r, theta, phi = sympy.symbols('r theta phi')
    zeta = sympy.symbols('zeta', positive=True, real=True)

    # --- Introduction and Problem Setup ---
    print("This script evaluates the integral <phi_i | 1/r | phi_j> for a 1s Slater orbital.")
    print("This means we are solving for <phi_1s| 1/r |phi_1s>.")
    print("The normalized 1s Slater-type orbital is given by:")
    print("phi_1s(r) = (zeta**3 / pi)**(1/2) * exp(-zeta*r)\n")
    print("The integral to solve in spherical coordinates is:")
    print("Integral[ (phi_1s)**2 * (1/r) * r**2 * sin(theta) ] dr d(theta) d(phi)\n")
    print("This integral can be broken down into four parts: a constant part from the")
    print("wavefunction's normalization, two angular integrals, and a radial integral.\n")

    # --- Step-by-step Calculation of Each Part ---
    print("--- Evaluating Each Part of the Integral ---")

    # Part 1: The constant part. This comes from squaring the normalization constant of phi_1s.
    # ( (zeta**3 / pi)**(1/2) )**2 = zeta**3 / pi
    constant_part = zeta**3 / sympy.pi
    print(f"1. Constant from normalization squared: {sympy.sstr(constant_part)}")

    # Part 2: The angular integral over the variable phi.
    # The integral of 1 from 0 to 2*pi.
    integral_phi = sympy.integrate(1, (phi, 0, 2 * sympy.pi))
    print(f"2. Angular integral over phi: Integral(d_phi) from 0 to 2*pi = {sympy.sstr(integral_phi)}")

    # Part 3: The angular integral over the variable theta.
    # The integral of sin(theta) from 0 to pi.
    integral_theta = sympy.integrate(sympy.sin(theta), (theta, 0, sympy.pi))
    print(f"3. Angular integral over theta: Integral(sin(theta)*d_theta) from 0 to pi = {sympy.sstr(integral_theta)}")

    # Part 4: The radial integral. The radial part of the integrand is found by combining
    # the radial parts of the wavefunctions, the operator, and the volume element:
    # (exp(-zeta*r))**2 * (1/r) * r**2 = r * exp(-2*zeta*r)
    radial_integrand = r * sympy.exp(-2 * zeta * r)
    integral_r = sympy.integrate(radial_integrand, (r, 0, sympy.oo))
    print(f"4. Radial integral: Integral(r * exp(-2*zeta*r)*dr) from 0 to infinity = {sympy.sstr(integral_r)}\n")

    # --- Assembling the Final Result ---
    print("--- Assembling the Final Equation ---")
    print("The value of the integral is the product of these four calculated parts.")

    # Convert each part to a string for clear printing
    part1_str = sympy.sstr(constant_part)
    part2_str = sympy.sstr(integral_phi)
    part3_str = sympy.sstr(integral_theta)
    part4_str = sympy.sstr(integral_r)

    # Calculate the final result by multiplying the parts
    final_result = constant_part * integral_phi * integral_theta * integral_r

    # Print the full equation with the value of each part
    print(f"Final Equation: ( {part1_str} ) * ( {part2_str} ) * ( {part3_str} ) * ( {part4_str} )")
    
    # Print the simplified final answer
    print(f"              = {final_result}")

if __name__ == '__main__':
    evaluate_slater_integral()