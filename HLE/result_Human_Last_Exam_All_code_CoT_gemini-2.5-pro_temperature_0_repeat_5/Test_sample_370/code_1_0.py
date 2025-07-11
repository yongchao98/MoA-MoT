import sympy as sp

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering
    mediated by a pseudoscalar boson in the low-energy limit (E << M).
    """
    # Step 1: Define symbolic variables
    # g: coupling constant
    # E: energy of one fermion in the CM frame
    # M: mass of the scalar mediator phi
    # theta, phi: scattering angles in the CM frame
    g, E, M, theta, phi = sp.symbols('g E M theta phi', real=True, positive=True)

    # Step 2: Explain the physical context and approximations
    print("This script calculates the total cross section for fermion-fermion scattering (psi + psi -> psi + psi).")
    print("The calculation is performed in the center-of-mass (CM) frame under the following approximations:")
    print("1. The fermion mass 'm' is negligible (high-energy fermions, E >> m).")
    print("2. The scattering energy is much smaller than the mediator mass 'M' (low-energy scattering, E << M).")
    print("-" * 50)

    # Step 3: Define Mandelstam variables in the CM frame for m=0
    # s: center-of-mass energy squared
    # t, u: momentum transfer squared
    s = 4 * E**2
    t = -2 * E**2 * (1 - sp.cos(theta))
    u = -2 * E**2 * (1 + sp.cos(theta))
    print("The Mandelstam variables (s, t, u) in the CM frame (for m=0) are:")
    print(f"s = {s}")
    print(f"t = {t}")
    print(f"u = {u}")
    print("-" * 50)

    # Step 4: Define the spin-averaged squared matrix element |M|^2
    # In the limit E << M, the propagators 1/(t-M^2) and 1/(u-M^2) are approx -1/M^2.
    # The full expression for |M|^2 simplifies to:
    M2_avg_approx = (g**4 / M**4) * (t**2 + u**2 - t * u)
    print("The spin-averaged squared matrix element |M|^2 in the E << M limit is:")
    print("|M|^2_avg = (g^4 / M^4) * (t^2 + u^2 - t*u)")
    
    # Simplify the expression for |M|^2
    M2_avg_simplified = sp.simplify(M2_avg_approx)
    print("\nSubstituting the expressions for t and u, this simplifies to:")
    print(f"|M|^2_avg = {M2_avg_simplified}")
    print("-" * 50)

    # Step 5: Define the differential cross section d(sigma)/d(Omega)
    # The formula is d(sigma)/d(Omega) = |M|^2_avg / (64 * pi^2 * s)
    dsig_dOmega = M2_avg_simplified / (64 * sp.pi**2 * s)
    dsig_dOmega_simplified = sp.simplify(dsig_dOmega)
    print("The differential cross section d(sigma)/d(Omega) is |M|^2_avg / (64 * pi^2 * s):")
    print(f"d(sigma)/d(Omega) = {dsig_dOmega_simplified}")
    print("-" * 50)

    # Step 6: Integrate over the solid angle to get the total cross section sigma
    # d(Omega) = sin(theta) d(theta) d(phi)
    # Integration limits: phi from 0 to 2*pi, theta from 0 to pi
    # We account for identical fermions in the final state by integrating theta from 0 to pi/2
    # and multiplying by 2pi for phi, then not dividing by the symmetry factor 2.
    # Or, integrate over the full range and divide by 2. Let's integrate over the full range for clarity.
    # The symmetry factor 1/2 for identical initial particles is already in the formula for sigma.
    integral_theta = sp.integrate(dsig_dOmega_simplified * sp.sin(theta), (theta, 0, sp.pi))
    total_sigma = sp.integrate(integral_theta, (phi, 0, 2 * sp.pi))
    
    print("To find the total cross section sigma, we integrate d(sigma)/d(Omega) over the full solid angle:")
    print("sigma = Integral( d(sigma)/d(Omega) * sin(theta) d(theta) d(phi) )")
    print("\nThe final result for the total cross section is:")
    
    # Step 7: Print the final equation with its components
    num, den = total_sigma.as_numer_denom()
    
    # Extract and print components for clarity as requested
    g_term = "g^4"
    E_term = "E^2"
    num_factor = "1" # The coefficient is 1
    
    den_factor = den.args[0]
    pi_term = "pi"
    M_term = "M^4"
    
    print(f"sigma = ({g_term} * {E_term}) / ({den_factor} * {pi_term} * {M_term})")

if __name__ == '__main__':
    calculate_cross_section()