import sympy as sp

def solve_integral():
    """
    Symbolically evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.
    """
    # Step 1: Define the necessary mathematical symbols.
    # r is the radial distance, zeta is the orbital exponent.
    # Both are positive real numbers.
    r, zeta = sp.symbols('r zeta', positive=True)

    # The 1s Slater orbital has the form N * exp(-zeta*r), where N is a normalization constant.
    # The square of the normalized orbital is (zeta**3/pi) * exp(-2*zeta*r).
    phi_1s_squared = (zeta**3 / sp.pi) * sp.exp(-2 * zeta * r)

    # Step 2: Define the operator.
    operator = 1/r

    # Step 3: Define the volume element in spherical coordinates.
    # dV = r**2 * sin(theta) * dr * dtheta * dphi.
    # The integral over the angular parts (dtheta, dphi) is 4*pi.
    angular_integral = 4 * sp.pi

    # Step 4: Set up the radial part of the integrand.
    # The full integrand is phi_1s_squared * operator * r**2.
    # We will integrate this with respect to r.
    radial_integrand = phi_1s_squared * operator * r**2

    # Step 5: Perform the definite integral for the radial part from 0 to infinity.
    radial_integral_result = sp.integrate(radial_integrand, (r, 0, sp.oo))

    # Step 6: The final result is the product of the radial integral and the angular integral.
    final_result = radial_integral_result * angular_integral

    # Step 7: Print the results in a step-by-step fashion.
    print("Evaluation of the integral <phi_1s| 1/r |phi_1s>")
    print("=" * 50)
    
    print(f"The squared normalized 1s Slater orbital is:\nphi_1s(r)^2 = {phi_1s_squared}\n")
    
    print(f"The operator is: O = {operator}\n")
    
    print("The integral is: Integral( phi_1s(r)^2 * O * dV ) over all space.")
    print("In spherical coordinates, dV = r^2 * sin(theta) * dr * dtheta * dphi.")
    print("The integral over the angular part is 4*pi.\n")

    print("The radial part of the integral is:")
    print(f"Integral from 0 to oo of ( {radial_integrand.simplify()} ) dr")
    
    # We calculate the result of the radial integral separately to show the steps.
    # The radial_integral_result is actually the integral of the radial_integrand / (4*pi)
    # The true radial integral is radial_integral_result * (4*pi) / (zeta**3/pi)
    # Let's show the integral of (r * exp(-2*zeta*r))
    radial_part = sp.integrate(r * sp.exp(-2*zeta*r), (r, 0, sp.oo))
    print(f"The core radial integral part evaluates to: 1 / (4 * zeta^2)")

    print("\nPutting it all together:")
    # The full expression is (zeta**3/pi) * (integral of r*exp(-2*zeta*r)) * (4*pi)
    # Which is (zeta**3/pi) * (1/(4*zeta**2)) * (4*pi)
    
    print("Final result = (Normalization_factor_squared) * (Radial_Integral) * (Angular_Integral)")
    print(f"Final result = ({zeta**3 / sp.pi}) * ({radial_part}) * ({angular_integral})")

    # Displaying the final equation with all its components
    print("\nThe final equation is:")
    print(f"⟨phi_1s| 1/r |phi_1s⟩ = ({sp.pretty(zeta**3 / sp.pi, use_unicode=False)}) * ({sp.pretty(radial_part, use_unicode=False)}) * (4*pi)")
    print(f"⟨phi_1s| 1/r |phi_1s⟩ = {final_result}")


if __name__ == '__main__':
    solve_integral()
