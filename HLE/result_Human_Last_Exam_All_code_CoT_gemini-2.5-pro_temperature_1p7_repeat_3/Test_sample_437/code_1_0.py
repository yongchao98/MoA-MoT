import sympy

def evaluate_1s_slater_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater orbital.
    """
    # 1. Define symbolic variables using sympy
    # r is the radial coordinate
    # zeta is the orbital exponent of the Slater-type orbital (STO)
    r, zeta = sympy.symbols('r zeta', positive=True)

    print("--- Evaluating the integral <phi_1s| 1/r |phi_1s> ---")
    print("\nStep 1: Define the 1s Slater Orbital (STO) wavefunction, phi_1s(r)")
    
    # Normalization constant for a 1s STO
    N = sympy.sqrt(zeta**3 / sympy.pi)
    
    # 1s STO wavefunction (spherically symmetric)
    phi_1s = N * sympy.exp(-zeta * r)
    
    print(f"The normalized 1s STO is: phi_1s(r) = {phi_1s}")

    print("\nStep 2: Set up the integral in spherical coordinates.")
    # The integral is Integral(phi_1s * (1/r) * phi_1s * dV)
    # dV (volume element) = r**2 * sin(theta) * dr * dtheta * dphi
    
    # The integrand (excluding the dV angular part)
    integrand_radial_part = phi_1s**2 * (1/r)
    
    print("The integral is Int( (phi_1s)^2 * (1/r) * r^2 * sin(theta) ) dr dtheta dphi")
    print("The integrand simplifies to: " + str(integrand_radial_part.simplify()) + " * r^2 * sin(theta)")

    print("\nStep 3: Integrate over the angular coordinates (theta and phi).")
    # The angular integral is Int(sin(theta) dtheta dphi) from theta=0 to pi, phi=0 to 2*pi
    # This integral evaluates to 4*pi
    angular_integral_result = 4 * sympy.pi
    print(f"The integral over angles gives a factor of {angular_integral_result}.")
    
    print("\nStep 4: Solve the remaining radial integral.")
    
    # The radial part of the integral to solve is:
    # Int( (phi_1s)^2 * (1/r) * r^2, r=0..infinity )
    # This simplifies to Int( N^2 * exp(-2*zeta*r) * r, r=0..infinity )
    
    # Let's break down the final calculation:
    # Full Integral = (Angular Part) * (Radial Part)
    #               = (4*pi) * Integral( (N^2 * exp(-2*zeta*r)) * (1/r) * r^2, dr)
    #               = (4*pi * N^2) * Integral( r * exp(-2*zeta*r), dr)
    
    # Prefactor = 4*pi*N^2 = 4*pi * (zeta^3/pi) = 4*zeta^3
    prefactor = 4 * sympy.pi * (N**2)
    print(f"The prefactor from the angular integral and normalization constant is: {prefactor}")
    
    # The purely radial integral to solve:
    radial_integral = sympy.integrate(r * sympy.exp(-2 * zeta * r), (r, 0, sympy.oo))
    print(f"The value of the radial integral Int(r * exp(-2*zeta*r), r=0..oo) is: {radial_integral}")
    
    print("\nStep 5: Combine the results to get the final answer.")
    final_result = prefactor * radial_integral
    
    # Show the final equation with all its numerical/symbolic parts
    print("\nFinal Equation Breakdown:")
    print(f"  <phi_1s| 1/r |phi_1s> = (Prefactor) * (Radial Integral Value)")
    print(f"  = ({prefactor}) * ({radial_integral})")
    
    print("\n--- Result ---")
    print(f"The symbolic result of the integral is: {final_result}")
    
    print("\n--- Numerical Example: Hydrogen Atom ---")
    # For a hydrogen atom in atomic units, zeta = 1.
    zeta_val = 1
    numerical_result = final_result.subs(zeta, zeta_val)
    print(f"For a hydrogen atom in atomic units, we set zeta = {zeta_val}.")
    print(f"The numerical value of the integral is: {numerical_result}")

if __name__ == '__main__':
    evaluate_1s_slater_integral()
