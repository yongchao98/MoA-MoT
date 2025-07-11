import math

def explain_npsd():
    """
    Explains and calculates the Noise Power Spectral Density (NPSD) for a simple reactor model.
    """
    # The double Fourier transform of the generalized pair correlation function
    # is commonly called the Noise Power Spectral Density (NPSD).
    function_name = "Noise Power Spectral Density (NPSD)"

    print(f"The space-time, double Fourier transform of the generalized pair correlation function is called the: {function_name}\n")

    # --- Model Parameters ---
    # We use a one-group diffusion model for a simple illustration.
    # The NPSD for this model has a Lorentzian shape determined by the
    # mode-dependent Rossi-alpha.

    # Reactivity in dollars (rho / beta_eff)
    rho_dollars = -0.5  # System is subcritical
    # Effective delayed neutron fraction
    beta_eff = 0.0065
    # Prompt neutron generation time (s)
    Lambda = 2.0e-5  # seconds
    # Neutron speed (cm/s)
    v_neutron = 2.2e5  # cm/s for thermal neutrons
    # Diffusion coefficient (cm)
    D = 0.5  # cm

    # Let's analyze a specific spatial mode 'k' and frequency 'omega'
    k = 0.1  # Spatial mode wave number (1/cm)
    omega = 300.0  # Angular frequency (rad/s)

    # --- Calculations ---
    # The reactivity in absolute units
    rho = rho_dollars * beta_eff

    # The fundamental mode Rossi-alpha (alpha for k=0)
    # alpha = (beta_eff - rho) / Lambda
    alpha_0 = (beta_eff - rho) / Lambda

    # The mode-dependent Rossi-alpha includes the spatial leakage term D*v*k^2
    # The term D*v*k^2 represents the increased decay rate of higher spatial modes.
    alpha_k = alpha_0 + D * v_neutron * (k**2)

    # The NPSD is proportional to 1 / (alpha_k^2 + omega^2)
    # The constant of proportionality depends on the source strength and detector efficiency.
    # We will represent it as 'C'.

    # --- Output ---
    print("For a one-group diffusion model, the NPSD has the form:")
    print("S(k, omega) = C / (alpha_k**2 + omega**2)")
    print("where alpha_k is the mode-dependent Rossi-alpha.\n")
    print("--- Example Calculation ---")
    print(f"Given parameters:")
    print(f"  Reactivity (rho): {rho:.6f} (or {rho_dollars} dollars)")
    print(f"  Effective delayed neutron fraction (beta): {beta_eff}")
    print(f"  Prompt neutron generation time (Lambda): {Lambda} s")
    print(f"  Diffusion coefficient (D): {D} cm")
    print(f"  Neutron speed (v): {v_neutron} cm/s\n")

    print(f"For the chosen spatial mode k = {k} 1/cm:")
    print(f"The mode-dependent Rossi-alpha (alpha_k) is calculated as:")
    print(f"alpha_k = (beta - rho)/Lambda + D*v*k^2")
    print(f"alpha_k = ({beta_eff:.4f} - ({rho:.6f}))/{Lambda:.1e} + {D}*{v_neutron:.1e}*({k})**2")
    print(f"alpha_k = {alpha_k:.2f} rad/s\n")

    print("Therefore, the final equation for the NPSD at the specified point is:")
    print(f"S(k={k}, omega={omega}) = C / ( ({alpha_k:.2f})**2 + ({omega:.2f})**2 )")


if __name__ == "__main__":
    explain_npsd()