import sympy
from sympy import Symbol, exp, pi, integrate, oo, sin, cos

def evaluate_1s_integral():
    """
    This function evaluates the integral <φ_1s| 1/r |φ_1s> for a Slater-type orbital.
    It breaks down the problem and shows the result of each part of the calculation.
    """
    # Define the symbols used in the calculation
    r = Symbol('r', positive=True)
    theta = Symbol('theta')
    phi = Symbol('phi')
    zeta = Symbol('zeta', positive=True)

    print("Evaluating the integral <φ_1s| 1/r |φ_1s> where φ_1s is a normalized 1s Slater orbital.")
    print("The integral in spherical coordinates is: I = ∫∫∫ (ζ³/π) * exp(-2ζr) * (1/r) * r² * sin(θ) dr dθ dφ")
    print("This simplifies to a product of three parts: I = (Constant) * (Radial Integral) * (Angular Integral)\n")

    # Part 1: The constant pre-factor
    # This comes from the square of the normalization constant of the 1s STO.
    pre_factor_str = "ζ³/π"
    pre_factor_sym = zeta**3 / pi
    print(f"1. The constant pre-factor from normalization is: {pre_factor_str}")

    # Part 2: The angular integral
    # The angular part of the integral is ∫ from 0 to 2π dφ * ∫ from 0 to π sin(θ) dθ
    angular_integral_phi = integrate(1, (phi, 0, 2*pi))
    angular_integral_theta = integrate(sin(theta), (theta, 0, pi))
    angular_integral_val = angular_integral_phi * angular_integral_theta
    print(f"2. The angular integral is (∫dφ) × (∫sin(θ)dθ) = ({angular_integral_phi}) × ({angular_integral_theta}) = {angular_integral_val}")

    # Part 3: The radial integral
    # The radial part of the integral is ∫ from 0 to ∞ of r * exp(-2ζr) dr
    # This comes from the original r-dependent parts: exp(-ζr) * (1/r) * exp(-ζr) * r² = r * exp(-2ζr)
    radial_integrand = r * exp(-2 * zeta * r)
    radial_integral_val = integrate(radial_integrand, (r, 0, oo))
    print(f"3. The radial integral is ∫ r·exp(-2ζr)dr from 0 to ∞, which evaluates to: {radial_integral_val}\n")

    # Combine all parts to form the final equation
    print("Combining all parts gives the final equation. Each term in the product is shown:")
    # Using sympy's pretty print for better visualization of the equation
    # The components are the "numbers" requested in the prompt.
    equation = sympy.Eq(sympy.Symbol("Value"), pre_factor_sym * angular_integral_val * radial_integral_val)
    print(f"   Value = ( {pre_factor_str} ) × ( {angular_integral_val} ) × ( {radial_integral_val} )")

    # Calculate the final result by simplifying the symbolic expression
    final_result = sympy.simplify(equation.rhs)

    print(f"\nSimplifying the expression leads to the final result:")
    print(final_result)


if __name__ == '__main__':
    evaluate_1s_integral()