import sympy as sp

def evaluate_slater_integral():
    """
    Evaluates the integral <phi_1s | 1/r | phi_1s> for a Slater-type orbital.
    """
    # 1. Define the symbols and the unnormalized 1s Slater orbital
    # r is the radial distance, zeta is the orbital exponent. Both are positive real numbers.
    r = sp.symbols('r', real=True, positive=True)
    zeta = sp.symbols('zeta', real=True, positive=True)

    # The unnormalized 1s Slater-type orbital (STO) is proportional to exp(-zeta*r)
    phi_unnormalized = sp.exp(-zeta * r)
    print(f"The unnormalized 1s Slater orbital is: phi_unnormalized = {phi_unnormalized}\n")

    # 2. Normalize the orbital
    # The normalization condition is Integral( (N*phi)^2 * d_tau ) = 1, where N is the normalization constant
    # and d_tau is the volume element in spherical coordinates (r^2 * sin(theta) dr dtheta dphi).
    # We can separate the integral into angular and radial parts.

    # The angular part of the integral is constant for any s-orbital.
    theta, phi = sp.symbols('theta phi', real=True)
    angular_integral = sp.integrate(sp.sin(theta), (theta, 0, sp.pi)) * sp.integrate(1, (phi, 0, 2*sp.pi))
    print(f"The angular part of the normalization integral evaluates to: {sp.sstr(angular_integral)}\n")

    # The radial part of the normalization integral for the unnormalized function:
    # Integral( phi_unnormalized**2 * r**2 dr) from 0 to infinity
    radial_integral_norm = sp.integrate(phi_unnormalized**2 * r**2, (r, 0, sp.oo))
    print(f"The radial part of the normalization integral is: {sp.sstr(radial_integral_norm)}\n")

    # The square of the normalization constant is N^2 = 1 / (angular_integral * radial_integral_norm)
    N_squared = 1 / (angular_integral * radial_integral_norm)
    N = sp.sqrt(N_squared)
    print(f"The square of the normalization constant is: N^2 = {sp.sstr(N_squared)}")
    phi_normalized = N * phi_unnormalized
    print(f"The normalized 1s Slater orbital is: phi_1s = {sp.sstr(phi_normalized)}\n")

    # 3. Evaluate the expectation value <phi_1s | 1/r | phi_1s>
    # The integral is Integral( phi_normalized**2 * (1/r) * d_tau )
    # This is equal to N^2 * Integral( phi_unnormalized**2 * (1/r) * d_tau )
    print("--- Evaluating the integral <phi_1s | 1/r | phi_1s> ---")
    
    # The angular part is the same as before.
    print(f"The angular part of the integral is again: {sp.sstr(angular_integral)}")

    # The radial part of the expectation value integral is:
    # Integral( phi_unnormalized**2 * (1/r) * r**2 dr ) = Integral( exp(-2*zeta*r) * r dr )
    radial_integral_exp = sp.integrate(phi_unnormalized**2 * (1/r) * r**2, (r, 0, sp.oo))
    print(f"The radial part of the expectation value integral is: {sp.sstr(radial_integral_exp)}\n")

    # 4. Combine the parts to get the final result
    # Result = N**2 * angular_integral * radial_integral_exp
    final_result = N_squared * angular_integral * radial_integral_exp

    print("The final result is the product of these three components:")
    print(f"<phi|1/r|phi> = (N^2) * (Angular Integral) * (Radial Integral)")
    
    # To satisfy the "output each number" requirement, we construct the equation string manually.
    N_squared_str = f"({sp.sstr(zeta)}**3/pi)"
    angular_int_str = f"(4*pi)"
    radial_int_exp_str = f"(1/(4*{sp.sstr(zeta)}**2))"
    print(f"The equation with values is: <phi|1/r|phi> = {N_squared_str} * {angular_int_str} * {radial_int_exp_str}")

    print("\nSimplifying the expression gives the final answer:")
    print(f"<phi_1s | 1/r | phi_1s> = {final_result}")

if __name__ == '__main__':
    evaluate_slater_integral()