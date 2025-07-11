import sympy

def evaluate_slater_integral():
    """
    Symbolically evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater-Type Orbital (STO).
    """
    # Define the necessary mathematical symbols with their properties
    # zeta is the orbital exponent, a positive real number.
    # r is the radial distance from the nucleus, a positive real number.
    zeta = sympy.Symbol('zeta', positive=True, real=True)
    r = sympy.Symbol('r', positive=True, real=True)
    pi = sympy.pi

    print("This script evaluates the expectation value of the 1/r operator for a 1s Slater-Type Orbital (STO).")
    print("The integral to be solved is: <phi_1s| 1/r |phi_1s>\n")

    # --- Step 1: Define the components of the integral ---
    print("Step 1: Define the components of the integral")
    # The square of the normalization constant for a 1s STO
    N_squared = zeta**3 / pi
    print(f"The square of the normalization constant, N^2, is: {N_squared}")

    # The angular part of the integral in spherical coordinates.
    # This comes from integrating sin(theta) over theta [0, pi] and dphi over phi [0, 2pi].
    # Integral[sin(theta)dtheta] from 0 to pi = 2
    # Integral[dphi] from 0 to 2*pi = 2*pi
    angular_part = 4 * pi
    print(f"The integral over the angular coordinates (theta, phi) evaluates to: {angular_part}")

    # The radial part of the integrand, after including r^2 from the volume element
    # and simplifying with the 1/r operator.
    # The expression inside the radial integral is: [exp(-zeta*r)]^2 * (1/r) * r^2 = r * exp(-2*zeta*r)
    radial_integrand = r * sympy.exp(-2 * zeta * r)
    print(f"The radial part of the integrand is: {radial_integrand}\n")


    # --- Step 2: Evaluate the radial integral ---
    print("Step 2: Evaluate the radial integral")
    # Integrate the radial part from r = 0 to infinity
    radial_integral_result = sympy.integrate(radial_integrand, (r, 0, sympy.oo))
    print(f"The value of the radial integral 'Integral({radial_integrand}) dr' from 0 to oo is: {radial_integral_result}\n")


    # --- Step 3: Combine all parts for the final result ---
    print("Step 3: Combine the components to get the final result")
    print("The total expectation value is the product of the three components: N^2 * (Angular Part) * (Radial Integral)")
    
    # Calculate the final result by multiplying the components
    final_result = N_squared * angular_part * radial_integral_result
    
    # Print the final equation with each component clearly shown
    print("\nFinal Equation:")
    print(f"   <phi_1s| 1/r |phi_1s> = ({N_squared}) * ({angular_part}) * ({radial_integral_result})")
    
    # Show simplification step
    simplified_term = N_squared * angular_part
    print(f"   Simplifying N^2 and the angular part: ({N_squared}) * ({angular_part}) = {simplified_term}")
    print(f"   So, the equation becomes: ({simplified_term}) * ({radial_integral_result})")
    
    # Print the final symbolic answer
    print(f"\n   The final result is: {final_result}")
    
    print("\n------------------------------------------------------------")
    print("Therefore, the value of the integral <phi_1s| 1/r |phi_1s> is simply the orbital exponent, zeta.")
    print("For a hydrogen atom in its ground state, zeta = 1, so the expectation value is 1.")
    print("------------------------------------------------------------")

if __name__ == '__main__':
    evaluate_slater_integral()
