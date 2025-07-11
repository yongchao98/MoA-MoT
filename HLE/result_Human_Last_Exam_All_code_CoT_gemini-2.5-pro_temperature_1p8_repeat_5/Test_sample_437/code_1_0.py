import sympy as sp

def evaluate_1s_integral():
    """
    Symbolically evaluates the integral <phi_1s| 1/r |phi_1s> for a 
    1s Slater-type orbital using the sympy library.

    It prints the components of the calculation and the final simplified result.
    """
    # Define the symbols needed for the calculation.
    # r, theta, phi are the spherical coordinates.
    # zeta is the orbital exponent of the Slater-type orbital.
    r, theta, phi = sp.symbols('r theta phi', real=True, positive=True)
    zeta = sp.symbols('zeta', real=True, positive=True)

    print("We are evaluating the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.")
    print("The integral is expressed as: Integral(phi_1s^* * (1/r) * phi_1s * dV)")
    print("where dV = r^2 * sin(theta) * dr * dtheta * dphi.")
    print("-" * 50)
    
    # The integral can be separated into the normalization constants and the integral itself.
    # Let phi = N * phi_unnormalized. The integral becomes N^2 * Integral(phi_u^* (1/r) phi_u dV).

    # First, let's determine the square of the normalization constant, N^2.
    # N^2 = 1 / Integral(phi_unnormalized^2 dV)
    # where phi_unnormalized = exp(-zeta*r)
    # The integrand for normalization is exp(-2*zeta*r) * r^2 * sin(theta)
    
    # The angular part of the normalization integral: Integral(sin(theta) dtheta dphi)
    angular_integral_norm = sp.integrate(sp.sin(theta), (phi, 0, 2 * sp.pi), (theta, 0, sp.pi))
    
    # The radial part of the normalization integral: Integral(exp(-2*zeta*r) * r^2 dr)
    radial_integral_norm = sp.integrate(sp.exp(-2 * zeta * r) * r**2, (r, 0, sp.oo))
    
    # Calculate N^2
    N_squared = 1 / (angular_integral_norm * radial_integral_norm)
    
    # Now, let's evaluate the main integral: Integral(exp(-zeta*r) * (1/r) * exp(-zeta*r) * dV)
    # This simplifies to: Integral(exp(-2*zeta*r) * r * sin(theta) * dr dtheta dphi)
    
    # The angular part of the main integral is the same as for normalization.
    angular_integral_main = angular_integral_norm
    
    # The radial part of the main integral: Integral(exp(-2*zeta*r) * r dr)
    radial_integral_main = sp.integrate(sp.exp(-2 * zeta * r) * r, (r, 0, sp.oo))

    # The final result is the product of N^2 and the integrals.
    final_result = N_squared * angular_integral_main * radial_integral_main
    
    # Print the "final equation" with each component, as requested.
    print("The final result is calculated from the following equation:")
    print("Result = (N^2) * (Angular Integral) * (Radial Integral)\n")
    print("The symbolic value of each component is:")
    print(f"  N^2 (Squared Normalization Constant) = {str(N_squared)}")
    print(f"  Angular Integral = {str(angular_integral_main)}")
    print(f"  Radial Integral for <1/r> = {str(radial_integral_main)}\n")

    print("Substituting these into the equation:")
    print(f"Result = ({str(N_squared)}) * ({str(angular_integral_main)}) * ({str(radial_integral_main)})\n")

    # The sp.simplify() function will perform the algebraic simplification.
    simplified_result = sp.simplify(final_result)
    
    print("After simplification, the final result of the integral is:")
    print(f"Result = {simplified_result}")


if __name__ == '__main__':
    # Execute the function to see the output.
    evaluate_1s_integral()