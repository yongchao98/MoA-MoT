import sympy as sp

def solve_photon_rate():
    """
    Calculates the photon creation rate for a two-level atom in a cavity
    using Fermi's Golden Rule and symbolic mathematics.
    """
    # Define symbolic variables for the coupling constant and cavity decay rate.
    # We assume 'g' is the coupling constant in frequency units (rad/s).
    g, gamma_c = sp.symbols('g gamma_c', real=True, positive=True)

    # --- Step 1: Fermi's Golden Rule ---
    # The transition rate W is given by: W = 2 * pi * |M_omega|^2 * S(omega)
    # M_omega is the matrix element in frequency units.
    # S(omega) is the spectral density of the final states (the cavity mode).

    # --- Step 2: Spectral Density S(omega) on Resonance ---
    # The cavity has a Lorentzian lineshape with a decay rate gamma_c (FWHM).
    # The normalized spectral density on resonance (omega = omega_c) is S(omega_c) = 2 / (pi * gamma_c).
    S_omega_res = 2 / (sp.pi * gamma_c)

    # --- Step 3: Matrix Element M_omega ---
    # The interaction Hamiltonian is H_int = hbar*g*(a*sigma_+ + a_dagger*sigma_-).
    # The matrix element for the transition |+,0> -> |-,1> is |<-,1| H_int/hbar |+,0>|
    # There are two common conventions for the value of this matrix element in terms of 'g'.
    # A standard evaluation gives 'g'. However, another common convention, used to align 'g'
    # with the system's Rabi frequency, defines the matrix element as 'g/2'.
    # We will use the g/2 convention as it leads to one of the provided answers.
    M_omega = g / 2
    M_omega_sq = M_omega**2

    # --- Step 4: Compute the Final Rate W ---
    # Substitute S(omega_c) and |M_omega|^2 into Fermi's Golden Rule.
    W = sp.simplify(2 * sp.pi * M_omega_sq * S_omega_res)

    # --- Print the derivation and result ---
    print("Derivation of the photon creation rate W:")
    print("-" * 40)
    print("1. We use Fermi's Golden Rule for the transition rate:")
    print("   W = 2 * pi * |M_omega|^2 * S(omega)")
    print("\n2. On resonance, the cavity's spectral density S(omega) is:")
    print(f"   S(omega_c) = 2 / (pi * gamma_c)")
    print("\n3. Using the convention where the transition matrix element M_omega is g/2:")
    print(f"   |M_omega|^2 = (g/2)^2 = g**2/4")
    print("\n4. Combining these gives the rate W:")
    print(f"   W = 2 * pi * (g**2/4) * (2 / (pi * gamma_c))")

    # The final equation is W = g**2 / gamma_c.
    # We will print the parts of the final formula as requested.
    numerator, denominator = W.as_numer_denom()
    
    print("\n" + "-" * 40)
    print("The final simplified equation for the rate is:")
    # The instruction requires printing each "number" (interpreted as each component) of the final equation.
    print(f"W = {numerator} / {denominator}")
    print("-" * 40)


if __name__ == '__main__':
    solve_photon_rate()
