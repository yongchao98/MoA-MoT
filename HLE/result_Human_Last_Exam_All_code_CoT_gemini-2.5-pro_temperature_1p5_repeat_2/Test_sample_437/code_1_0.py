import sympy
from sympy import symbols, integrate, pi, oo, exp, sin, sqrt, simplify, Eq, Integral, pretty

def evaluate_1s_integral():
    """
    This function symbolically evaluates the integral <phi_1s| 1/r |phi_1s>
    for a Slater-type 1s orbital and prints the derivation steps.
    """
    # Define the symbols to be used in the calculation.
    # r: radial distance
    # theta: polar angle
    # phi: azimuthal angle
    # zeta: orbital exponent (a positive real number)
    r, theta, phi = symbols('r theta phi', real=True)
    zeta = symbols('zeta', positive=True)

    # The integral I = <phi_1s| 1/r |phi_1s> can be broken down into three parts:
    # I = (Normalization Constant^2) * (Angular Integral) * (Radial Integral)
    # The script will calculate each part and then combine them.

    print("The final equation is of the form: I = (N^2) * (Angular Integral) * (Radial Integral)")
    print("-" * 70)

    # Part 1: The square of the normalization constant for a 1s STO.
    # The normalized 1s STO is phi_1s = N * exp(-zeta*r), where N = (zeta**3 / pi)**(1/2).
    norm_const_sq = zeta**3 / pi
    print("Part 1: The Square of the Normalization Constant (N^2)")
    print(f"The normalization constant N is sqrt(zeta^3 / pi).")
    print(f"N^2 = (sqrt(zeta^3 / pi))^2 = {pretty(norm_const_sq)}")
    print("-" * 70)

    # Part 2: The angular integral.
    # The operator and wavefunction are spherically symmetric, so the angular part of
    # the integral is over the volume element sin(theta).
    # Integral from 0 to 2*pi dphi, and integral from 0 to pi of sin(theta) dtheta.
    integral_over_phi = integrate(1, (phi, 0, 2*pi))
    integral_over_theta = integrate(sin(theta), (theta, 0, pi))
    angular_integral_result = integral_over_phi * integral_over_theta
    print("Part 2: The Angular Integral")
    print("The angular part is Integral( dphi ) * Integral( sin(theta) dtheta )")
    print(f"Integral over phi from 0 to 2*pi of 1 dphi = {integral_over_phi}")
    print(f"Integral over theta from 0 to pi of sin(theta) dtheta = {integral_over_theta}")
    print(f"Total angular integral result = {integral_over_phi} * {integral_over_theta} = {angular_integral_result}")
    print("-" * 70)

    # Part 3: The radial integral.
    # The radial part of the integral is Integral( |psi_unnormalized|^2 * (1/r) * r^2 dr )
    # psi_unnormalized = exp(-zeta*r)
    # So the integrand is (exp(-2*zeta*r)) * (1/r) * r^2 = r * exp(-2*zeta*r).
    radial_integrand = r * exp(-2 * zeta * r)
    radial_integral_result = integrate(radial_integrand, (r, 0, oo))
    print("Part 3: The Radial Integral")
    radial_integral_expr = Integral(radial_integrand, (r, 0, oo))
    print(f"The radial integral is: {pretty(radial_integral_expr)}")
    print(f"Result of the radial integral = {pretty(radial_integral_result)}")
    print("-" * 70)

    # Combine all parts to get the final result.
    final_result = norm_const_sq * angular_integral_result * radial_integral_result
    
    print("Final Calculation: Combining all parts")
    print("I = N^2 * (Angular Integral) * (Radial Integral)")
    # We print each numerical or symbolic part of the final equation.
    print(f"I = ({pretty(norm_const_sq)}) * ({pretty(angular_integral_result)}) * ({pretty(radial_integral_result)})")
    
    # Let sympy simplify the final expression.
    simplified_result = simplify(final_result)
    print(f"\nAfter simplification, the evaluation of the integral <phi_1s| 1/r |phi_1s> is:")
    print(f"I = {simplified_result}")


# Execute the function to perform the calculation and print the results.
evaluate_1s_integral()