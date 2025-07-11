import sympy as sp

def evaluate_slater_integral():
    """
    This function symbolically evaluates the integral <phi | 1/r | phi> for a 
    normalized 1s Slater-type orbital (STO). This integral corresponds to the
    expectation value of the electron-nucleus potential energy.
    """
    
    # Step 1: Define the necessary mathematical symbols.
    # r is the radial distance from the nucleus.
    # zeta is the orbital exponent, a positive constant that determines the orbital's size.
    # theta and phi are the standard angular coordinates in a spherical system.
    r, zeta = sp.symbols('r zeta', positive=True)
    theta, phi = sp.symbols('theta phi', real=True)

    # Step 2: Define the normalized 1s Slater orbital wavefunction.
    # The function is given by phi_1s = N * exp(-zeta*r), where N is a normalization constant.
    # For a 1s STO, the normalization constant N is (zeta^3 / pi)^(1/2).
    # This ensures that the probability of finding the electron anywhere in space is 1.
    # Integral(|phi_1s|^2 dV) = 1
    normalization_constant = sp.sqrt(zeta**3 / sp.pi)
    phi_1s = normalization_constant * sp.exp(-zeta * r)

    # Step 3: Set up the integral for the expectation value.
    # The integral is <phi_1s | 1/r | phi_1s> = Integral(phi_1s^* * (1/r) * phi_1s dV)
    # Since phi_1s is real, the complex conjugate phi_1s^* is just phi_1s.
    operator = 1 / r
    
    # The volume element in spherical coordinates is dV = r^2 * sin(theta) dr dtheta dphi.
    volume_element = r**2 * sp.sin(theta)

    # The complete expression to be integrated is:
    integrand = phi_1s * operator * phi_1s * volume_element

    # Step 4: Evaluate the definite triple integral.
    # SymPy can compute this integral over all space:
    # r from 0 to infinity
    # theta from 0 to pi
    # phi from 0 to 2*pi
    integral_result = sp.integrate(integrand, (r, 0, sp.oo), (theta, 0, sp.pi), (phi, 0, 2*sp.pi))

    # Step 5: Display the final result.
    # The result of the evaluation is the final term in the equation.
    lhs = "<phi_1s | 1/r | phi_1s>"
    
    print("The integral to be evaluated is <phi | 1/r | phi> for a 1s Slater orbital.")
    print(f"The normalized 1s Slater orbital is defined as: phi_1s = {normalization_constant} * exp(-zeta*r)")
    print(f"The operator is: {operator}")
    print("\nThe final equation is:")
    # The instruction asked to show numbers in the final equation. 
    # Here the result is symbolic. The equation itself is the output.
    print(f"{lhs} = {integral_result}")

if __name__ == '__main__':
    evaluate_slater_integral()