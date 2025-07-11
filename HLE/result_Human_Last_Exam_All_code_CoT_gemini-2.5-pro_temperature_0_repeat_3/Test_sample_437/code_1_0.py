import sympy
from sympy import symbols, exp, integrate, pi, oo, sin, pretty_print, Eq, Symbol

def evaluate_slater_integral():
    """
    Symbolically evaluates the integral <phi_1s | 1/r | phi_1s> for a 
    Slater-type orbital using the sympy library.
    """
    # Define the symbolic variables. zeta is the orbital exponent.
    r, theta, phi = symbols('r theta phi', real=True)
    zeta = symbols('zeta', positive=True, real=True)
    
    print("This script evaluates the integral <phi_1s | 1/r | phi_1s>.\n")
    print("The integral is separated into three parts:")
    print("Result = (Normalization_Factor) * (Radial_Integral) * (Angular_Integral)\n")

    # Part 1: The normalization factor
    # The normalized 1s STO is phi_1s = sqrt(zeta**3 / pi) * exp(-zeta*r).
    # The factor from the squared wavefunctions is (phi_1s)^2 / (exp(-zeta*r))^2
    norm_factor = zeta**3 / pi
    print("1. The normalization factor from the squared wavefunctions is:")
    pretty_print(Eq(Symbol('Normalization_Factor'), norm_factor))
    print("-" * 40)

    # Part 2: The angular integral
    # The angular part of the integral is over d(phi) and sin(theta)d(theta)
    angular_integral_phi = integrate(1, (phi, 0, 2*pi))
    angular_integral_theta = integrate(sin(theta), (theta, 0, pi))
    angular_integral = angular_integral_phi * angular_integral_theta
    
    print("2. The angular integral is Integral(sin(theta) d(theta) d(phi))")
    print("   (from theta=0 to pi, phi=0 to 2*pi). The result is:")
    pretty_print(Eq(Symbol('Angular_Integral'), angular_integral))
    print("-" * 40)

    # Part 3: The radial integral
    # The radial part is Integral(r * exp(-2*zeta*r) dr)
    # The 'r' comes from (1/r operator) * (r^2 from volume element)
    radial_integrand = r * exp(-2 * zeta * r)
    radial_integral = integrate(radial_integrand, (r, 0, oo))
    
    print("3. The radial integral is Integral(r * exp(-2*zeta*r) dr)")
    print("   (from r=0 to infinity). The result is:")
    pretty_print(Eq(Symbol('Radial_Integral'), radial_integral))
    print("-" * 40)

    # Final Calculation
    # Multiply the three parts together to get the final result.
    final_result = norm_factor * radial_integral * angular_integral
    
    print("Finally, we multiply these three parts together.")
    print("The final equation is:")
    # We use f-string formatting to display the components before simplification
    print(f"Result = ({norm_factor}) * ({radial_integral}) * ({angular_integral})")
    
    print("\nWhich simplifies to the final answer:")
    pretty_print(Eq(Symbol('Result'), final_result))

if __name__ == '__main__':
    evaluate_slater_integral()