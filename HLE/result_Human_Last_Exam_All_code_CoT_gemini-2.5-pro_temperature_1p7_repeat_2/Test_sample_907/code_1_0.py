import numpy as np
from scipy.constants import hbar, e

def get_absorption_equations():
    """
    This script calculates and prints the equations for the absorption cross-section
    for a molecular chain interacting with a Gaussian laser pulse, for two cases:
    a) non-interacting molecules
    b) molecules with near-neighbor interactions.
    """

    # --- Define Model Parameters ---
    # Number of molecules in the chain
    N = 10
    # Single-molecule excitation energy (in electron-volts)
    E_site_eV = 2.0
    # Inter-molecular coupling constant (in electron-volts).
    # A negative value is typical for J-aggregates.
    J_eV = -0.1
    # Transition dipole moment of a single molecule (in Debye). Let's assume its magnitude.
    # The orientation is handled by a constant factor we will omit.
    mu_eg_abs_debye = 5.0
    # Gaussian laser pulse duration (FWHM duration is approx. 1.665 * T)
    T_fs = 50.0

    # --- Convert parameters to SI units ---
    E_site_J = E_site_eV * e
    J_J = J_eV * e
    # We are mainly interested in the shape, so we can use |mu|^2 as a prefactor.
    # For simplicity, we'll set the numerical value of |mu|^2 to 1 in the final
    # expression and note that it scales the overall intensity.
    mu_sq_factor = mu_eg_abs_debye**2
    # The term we need for the exponent is T^2/hbar^2
    T_s = T_fs * 1e-15
    exp_factor = T_s**2 / (2 * hbar**2) # This is 1/(2*sigma_E^2) where sigma_E is energy uncertainty
    
    print("This program provides the equations for the absorption cross-section \u03C3(\u03C9_L)")
    print("of a molecular chain interacting with a Gaussian laser pulse with central frequency \u03C9_L.")
    print("The cross-section is proportional to the expression shown.")
    print("\n--- Based on the following parameters ---")
    print(f"Number of molecules N = {N}")
    print(f"Single molecule transition energy E_site = {E_site_eV} eV")
    print(f"Near-neighbor coupling J = {J_eV} eV")
    print(f"Pulse duration parameter T = {T_fs} fs")
    print(f"The symbol \u03C9_L represents the variable laser frequency in rad/s.")
    print("-" * 40)


    # --- Case a) No interaction (J=0) ---
    print("\n\na) Equation for non-interacting molecules (J = 0):")
    print("The molecules are independent. The absorption spectrum consists of a single peak")
    print("at the molecular transition frequency. The peak has a Gaussian shape due to the")
    print("finite duration of the laser pulse.")
    
    # Calculate transition frequency
    omega_fi = E_site_J / hbar
    # Calculate the exponent term T^2/2
    exp_factor_a = T_s**2 / 2
    
    print("\nFinal Equation (a):")
    # Output each number in the final equation!
    print(f"\u03C3_a(\u03C9_L) \u221D {mu_sq_factor:.2f} * exp( -({omega_fi:.3e} - \u03C9_L)\u00B2 * {exp_factor_a:.3e} )")
    print("\nWhere:")
    print(f"  - The prefactor {mu_sq_factor:.2f} represents the squared transition dipole moment magnitude (in Debye^2).")
    print(f"  - The transition frequency is \u03C9_fi = {omega_fi:.3e} rad/s (equivalent to {E_site_eV} eV).")
    print(f"  - The Gaussian width is determined by the pulse duration parameter T.")
    print(f"  - The expression in the exponent is -(\u03C9_fi - \u03C9_L)\u00B2 * T\u00B2/2.")

    print("\n" + "="*70 + "\n")

    # --- Case b) Near-neighbor interaction ---
    print("b) Equation for interacting molecules (J \u2260 0):")
    print("Interactions create delocalized exciton states with a band of energies.")
    print("Absorption occurs to specific states dictated by a selection rule (only k=odd states).")
    print("The total cross-section is a sum of Gaussian peaks, each centered at an allowed exciton transition frequency.")

    print("\nFinal Equation (b):")
    print(f"\u03C3_b(\u03C9_L) \u221D \u03A3 [ Oscillator_Strength_k * exp( -(\u03C9_k - \u03C9_L)\u00B2 * {exp_factor_a:.3e} ) ]")
    print("\nThe sum is over the allowed exciton states k. The terms are:")

    for k in range(1, N + 1):
        # Selection rule: only k=odd states are optically allowed for a simple linear chain
        if k % 2 != 0:
            # Calculate exciton energy E_k and frequency omega_k
            E_k_J = E_site_J + 2 * J_J * np.cos(np.pi * k / (N + 1))
            E_k_eV = E_k_J / e
            omega_k = E_k_J / hbar

            # Calculate the oscillator strength factor (proportional to |<k|mu|g>|^2)
            # The sum is |sum_n sin(pi*k*n/(N+1))|^2
            x = np.pi * k / (N + 1)
            # using analytical formula for sum of sines: sin(Nx/2)sin((N+1)x/2) / sin(x/2)
            # sin((N+1)x/2) = sin(k*pi/2), which is +/-1 for odd k
            sum_sin_term = np.sin(N*x/2) * np.sin((N+1)*x/2) / np.sin(x/2)
            # Oscillator strength is proportional to sum_sin_term^2
            # Let's normalize by the strongest peak (k=1)
            if k==1:
                norm_factor = sum_sin_term**2
            osc_strength = mu_sq_factor * (sum_sin_term**2 / norm_factor)
            
            # Output each number in the final equation!
            print(f"\n--- Term for k = {k} ---")
            print(f"  Exciton Energy E_{k} = {E_k_eV:.4f} eV")
            print(f"  Transition Frequency \u03C9_{k} = {omega_k:.3e} rad/s")
            print(f"  Term_{k}: {osc_strength:.3f} * exp( -({omega_k:.3e} - \u03C9_L)\u00B2 * {exp_factor_a:.3e} )")

if __name__ == '__main__':
    get_absorption_equations()
