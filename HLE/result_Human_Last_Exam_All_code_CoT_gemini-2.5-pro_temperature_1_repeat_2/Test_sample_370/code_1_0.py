import sympy as sp

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering
    mediated by a heavy scalar particle.
    """
    # Define symbolic variables
    g, E, M, theta, phi = sp.symbols('g E M theta phi', real=True, positive=True)
    pi = sp.pi

    print("Step 1: Define Kinematics and Approximations")
    print("---------------------------------------------")
    print("Process: Two-fermion scattering psi(p1) + psi(p2) -> psi(p3) + psi(p4)")
    print("Assumptions:")
    print("1. Fermions are massless (m=0), which is the high-energy limit (E >> m).")
    print("2. The scalar mediator is very heavy (M >> E), so the interaction is point-like.\n")

    # In the center-of-mass (CM) frame with m=0:
    s = 4 * E**2
    t_expr = -2 * E**2 * (1 - sp.cos(theta))
    u_expr = -2 * E**2 * (1 + sp.cos(theta))

    print("Step 2: Calculate Spin-Averaged Squared Amplitude <|M|^2>")
    print("---------------------------------------------------------")
    print("The squared amplitude is the sum of t-channel, u-channel, and interference terms.")
    print("In the heavy mediator limit, this simplifies significantly.")

    # In the limit M >> E, t, u, the propagators 1/(t-M^2) and 1/(u-M^2) both become approx -1/M^2.
    # The spin-averaged squared amplitude becomes:
    # <|M|^2> ~= (g^4 / M^4) * (t^2 + u^2 + t*u)
    
    t, u = sp.symbols('t u')
    M2_avg_general = (g**4 / M**4) * (t**2 + u**2 + t*u)
    print(f"Approximated <|M|^2> = {M2_avg_general}")

    # Substitute the CM frame expressions for t and u
    M2_avg_cm = M2_avg_general.subs({t: t_expr, u: u_expr})
    M2_avg_cm_simplified = sp.simplify(M2_avg_cm)

    print("\nSubstituting CM frame kinematics (s, t, u):")
    print(f"<|M|^2> = {M2_avg_cm_simplified}\n")

    print("Step 3: Calculate the Differential Cross Section d(sigma)/d(Omega)")
    print("----------------------------------------------------------------")
    # The general formula in the CM frame is d(sigma)/d(Omega) = (1 / (128 * pi^2 * s)) * <|M|^2>
    # The factor 128 = 2 * 64, where the extra 2 in the denominator compared to some textbooks
    # accounts for the flux factor for identical particles. When integrating to get the total
    # cross section, we must still account for identical final states.
    dsig_dOmega = (1 / (128 * pi**2 * s)) * M2_avg_cm_simplified
    dsig_dOmega_simplified = sp.simplify(dsig_dOmega)

    print("Using the formula d(sigma)/d(Omega) = (1 / (128 * pi^2 * s)) * <|M|^2>")
    print(f"d(sigma)/d(Omega) = {dsig_dOmega_simplified}\n")

    print("Step 4: Calculate the Total Cross Section sigma")
    print("---------------------------------------------")
    print("To get the total cross section for identical final particles, we integrate over half the solid angle,")
    print("or equivalently, integrate over the full solid angle and multiply by 1/2.")
    print("sigma = (1/2) * Integral( d(sigma)/d(Omega) * sin(theta) d(theta) d(phi) )")

    # Integrate over d(Omega) = sin(theta) d(theta) d(phi)
    # phi from 0 to 2*pi, theta from 0 to pi
    integrand = dsig_dOmega_simplified * sp.sin(theta)
    
    # Perform integration
    # The factor of 1/2 is for identical particles in the final state.
    sigma = sp.Rational(1, 2) * sp.integrate(integrand, (phi, 0, 2*pi), (theta, 0, pi))
    sigma_simplified = sp.simplify(sigma)
    
    print("\nFinal Result:")
    print("------------")
    # To satisfy the output format requirement "output each number in the final equation!"
    # we manually construct the string from the simplified sympy expression.
    num, den = sigma_simplified.as_numer_denom()
    
    # Extract numerical coefficients
    num_coeff = num.as_coeff_mul()[0]
    den_coeff = den.as_coeff_mul()[0]
    
    # Reconstruct string with explicit numbers
    num_str = str(num).replace(str(num_coeff), str(int(num_coeff)))
    den_str = str(den).replace(str(den_coeff), str(int(den_coeff)))
    
    print(f"The total cross section sigma is:")
    print(f"sigma = ({num_str}) / ({den_str})")
    print("\nExplicitly showing the numbers in the equation:")
    print(f"sigma = ( {int(num_coeff)} * g**4 * E**2 ) / ( {int(den_coeff)} * pi * M**4 )")


if __name__ == '__main__':
    calculate_cross_section()