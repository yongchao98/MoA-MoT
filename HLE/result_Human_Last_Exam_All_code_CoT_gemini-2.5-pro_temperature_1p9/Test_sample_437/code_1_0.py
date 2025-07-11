import sympy
from sympy import symbols, exp, pi, integrate, sin, oo, sqrt, pretty_print

def evaluate_1s_integral():
    """
    This function evaluates the integral <phi_1s| 1/r |phi_1s> for a 
    normalized 1s Slater-Type Orbital (STO) using symbolic mathematics.
    """
    print("Starting evaluation of the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.")
    print("-" * 70)
    
    # Define symbolic variables
    # r is the radial coordinate
    # zeta is the orbital exponent
    # theta and phi are the angular coordinates
    r, zeta = symbols('r zeta', positive=True)
    theta, phi = symbols('theta phi')

    # The normalized 1s Slater orbital is phi_1s = (zeta**3 / pi)**(1/2) * exp(-zeta*r)
    # The normalization constant squared is (zeta**3 / pi)
    norm_const_sq = zeta**3 / pi
    
    # The integral can be separated into three parts:
    # I = (Norm_Const^2) * (Radial Part) * (Angular Part)

    # 1. Normalization constant part
    # This is the square of the normalization constant of the 1s STO.
    part1_norm = norm_const_sq
    
    # 2. Radial part of the integral
    # The radial integrand is r * exp(-2*zeta*r)
    # This comes from (e^{-zeta*r}) * (1/r) * (e^{-zeta*r}) * r^2
    radial_integrand = r * exp(-2 * zeta * r)
    part2_radial = integrate(radial_integrand, (r, 0, oo))

    # 3. Angular part of the integral
    # The angular integrand is sin(theta)
    angular_integrand = sin(theta)
    # Integral over dphi from 0 to 2*pi gives 2*pi
    # Integral over sin(theta) dtheta from 0 to pi gives 2
    part3_angular = integrate(angular_integrand, (theta, 0, pi), (phi, 0, 2*pi))

    # The final result is the product of the three parts
    final_result = part1_norm * part2_radial * part3_angular
    
    print("The integral is a product of three terms:")
    print("\n1. Normalization Constant Squared:")
    pretty_print(part1_norm)
    
    print("\n2. Radial Integral: Integral from 0 to oo of (r * exp(-2*zeta*r)) dr")
    pretty_print(part2_radial)
    
    print("\n3. Angular Integral: Integral over all solid angles of (sin(theta)) d(theta) d(phi)")
    pretty_print(part3_angular)
    
    print("\n" + "-" * 70)
    print("The final equation is the product of these three parts:")
    # We print each "number" (symbolic part) of the final equation
    print(f"({part1_norm}) * ({part2_radial}) * ({part3_angular}) = {final_result}")
    print("\n" + "-" * 70)

    print("\nTherefore, the evaluated integral <phi_1s| 1/r |phi_1s> is:")
    pretty_print(final_result)

if __name__ == '__main__':
    evaluate_1s_integral()
