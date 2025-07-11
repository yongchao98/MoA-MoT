import sympy
from sympy import symbols, exp, pi, integrate, oo

def evaluate_slater_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital
    using symbolic mathematics.
    """
    # Define symbolic variables for the calculation
    r, zeta = symbols('r zeta', real=True, positive=True)

    # The integral is the product of three terms:
    # 1. The square of the normalization constant (N^2)
    # 2. The angular integral part
    # 3. The radial integral part

    # Part 1: The square of the normalization constant for a 1s STO.
    # N = (zeta**3 / pi)**(1/2), so N**2 is zeta**3 / pi.
    n_squared = zeta**3 / pi

    # Part 2: The integral over the angular coordinates (theta and phi).
    # Integral(sin(theta) d(theta) d(phi)) over the whole sphere is 4*pi.
    angular_integral = 4 * pi

    # Part 3: The radial part of the integral.
    # This comes from Integral( [exp(-zeta*r)] * (1/r) * [exp(-zeta*r)] * r**2 dr )
    # which simplifies to Integral( r * exp(-2*zeta*r) dr ) from 0 to infinity.
    radial_integrand = r * exp(-2 * zeta * r)
    radial_integral = integrate(radial_integrand, (r, 0, oo))

    # The final result is the product of these three components.
    final_result = n_squared * radial_integral * angular_integral

    # Print the breakdown of the calculation as requested.
    print("Evaluating the integral <phi_1s| 1/r |phi_1s> = (N²) * (Radial Integral) * (Angular Integral)")
    print("-" * 70)
    print("The values of the components are:")
    print(f"  - N² (Square of Normalization Constant): {n_squared}")
    print(f"  - Radial Integral: {radial_integral}")
    print(f"  - Angular Integral: {angular_integral}")
    print("-" * 70)

    # Print the final equation with all the numbers/symbols.
    term1_str = f"({str(n_squared)})"
    term2_str = f"({str(radial_integral)})"
    term3_str = f"({str(angular_integral)})"
    result_str = str(final_result)

    print("Final Equation:")
    print(f"{term1_str} * {term2_str} * {term3_str} = {result_str}")


if __name__ == '__main__':
    evaluate_slater_integral()