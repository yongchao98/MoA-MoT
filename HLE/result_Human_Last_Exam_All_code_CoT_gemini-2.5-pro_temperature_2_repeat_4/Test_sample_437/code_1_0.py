import sympy as sp

def evaluate_1s_slater_integral():
    """
    This function evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater-type orbital
    using the sympy library for symbolic mathematics.
    """
    # Define the symbols used in the integral
    # r: radial distance
    # theta: polar angle
    # phi: azimuthal angle
    # zeta: orbital exponent (a positive real number)
    r, theta, phi = sp.symbols('r theta phi', real=True)
    zeta = sp.symbols('zeta', positive=True)

    # The problem is to evaluate <phi_i | 1/r | phi_j> for a 1s Slater orbital.
    # This means phi_i = phi_j = phi_1s.
    # The normalized 1s Slater orbital is phi_1s = (zeta^3 / pi)^(1/2) * exp(-zeta*r)
    phi_1s_normalization_constant = sp.sqrt(zeta**3 / sp.pi)
    phi_1s = phi_1s_normalization_constant * sp.exp(-zeta * r)

    # The operator is 1/r
    operator = 1 / r

    # The integrand is phi_1s^* * operator * phi_1s * dV
    # where dV is the volume element in spherical coordinates: r^2 * sin(theta) dr d(theta) d(phi)
    # Since phi_1s is real, phi_1s^* = phi_1s
    integrand = phi_1s * operator * phi_1s * r**2 * sp.sin(theta)

    # Simplify the integrand before integration
    integrand_simplified = sp.simplify(integrand)
    
    # Perform the triple integration over all space:
    # r from 0 to infinity
    # theta from 0 to pi
    # phi from 0 to 2*pi
    integral_value = sp.integrate(
        integrand_simplified,
        (r, 0, sp.oo),
        (theta, 0, sp.pi),
        (phi, 0, 2*sp.pi)
    )

    # Output the result
    # We construct the equation string to be clear about what was calculated.
    # The components are the orbital phi_1s, the operator 1/r, and the resulting value.
    phi_i_str = "phi_1s"
    operator_str = "1/r"
    phi_j_str = "phi_1s"

    print(f"Evaluating the integral: <{phi_i_str} | {operator_str} | {phi_j_str}>")
    print(f"Where phi_1s(r) = (zeta^3 / pi)^(1/2) * exp(-zeta*r)")
    print("-" * 30)
    # In the final equation, we explicitly show each component.
    print(f"Final Equation: <{phi_i_str} | {operator_str} | {phi_j_str}> = {integral_value}")

if __name__ == '__main__':
    evaluate_1s_slater_integral()