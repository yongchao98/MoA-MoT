import numpy as np

def calculate_lifetime_and_ratio():
    """
    Calculates the theoretical lifetime of the Sodium 3p state and its ratio
    to the experimental value based on the provided hydrogenic model.
    """
    # --- Step 1: Define Constants and Parameters ---
    # Physical constants
    e = 1.602e-19  # Electron charge (C)
    a0 = 5.292e-11  # Bohr radius (m)
    c = 2.998e8     # Speed of light (m/s)
    h_bar = 1.0546e-34 # Reduced Planck constant (J*s)
    eps0 = 8.854e-12   # Vacuum permittivity (F/m)
    
    # Given experimental parameters
    lambda_exp = 589e-9  # Wavelength (m)
    tau_exp = 16.2e-9    # Experimental lifetime (s)
    
    # Model parameters
    # The problem provides Z=11 but asks to use hydrogenic wavefunctions.
    # For a realistic model of Sodium's outer electron, the effective
    # nuclear charge Z_eff is approximately 1 due to screening by inner electrons.
    Z_eff = 1.0
    n = 3  # Principal quantum number
    L_upper = 1 # Orbital quantum number for p-state
    S = 0.5 # Electron spin

    print("--- Input Parameters ---")
    print(f"Wavelength (λ): {lambda_exp:.3e} m")
    print(f"Experimental Lifetime (τ_exp): {tau_exp:.3e} s")
    print(f"Effective Nuclear Charge (Z_eff): {Z_eff}")
    print("-" * 20 + "\n")

    # --- Step 2: Calculate Intermediate Quantities ---
    # Angular frequency of the transition
    omega = 2 * np.pi * c / lambda_exp
    
    # Radial integral for a hydrogenic n,L -> n,L-1 transition
    I_r = ( (3 * n) / (2 * Z_eff) ) * np.sqrt(n**2 - L_upper**2) * a0

    # Line strength S_21 = e^2 * (2S+1) * max(L1, L2) * I_r^2
    # max(L_upper, L_lower) = max(1, 0) = 1
    line_strength = (e**2) * (2*S + 1) * L_upper * (I_r**2)
    
    # Degeneracy of the upper state (3p) including spin
    g2 = (2 * L_upper + 1) * (2 * S + 1)
    
    print("--- Intermediate Calculations ---")
    print(f"Angular Frequency (ω = 2πc/λ): {omega:.4e} rad/s")
    print(f"Radial Integral (I_r): {I_r:.4e} m")
    print(f"Upper State Degeneracy (g2): {g2}")
    print(f"Line Strength (S): {line_strength:.4e} C^2 m^2")
    print("-" * 20 + "\n")
    
    # --- Step 3: Calculate Theoretical Lifetime ---
    # Einstein A coefficient (transition rate)
    # A = (ω^3 * S) / (3 * π * ε_0 * ħ * c^3 * g_2)
    numerator_A = omega**3 * line_strength
    denominator_A = 3 * np.pi * eps0 * h_bar * c**3 * g2
    A_coeff = numerator_A / denominator_A
    
    # Theoretical lifetime is the reciprocal of the A coefficient
    tau_th = 1 / A_coeff
    
    print("--- Lifetime Calculation ---")
    print(f"Einstein A Coefficient (A = ω³S / (3πε₀ħc³g₂))")
    print(f"A = ({omega:.3e}³ * {line_strength:.3e}) / (3π * {eps0:.3e} * {h_bar:.3e} * {c:.3e}³ * {g2})")
    print(f"A = {A_coeff:.4e} s⁻¹")
    print(f"Theoretical Lifetime (τ_th = 1/A): {tau_th:.4e} s")
    print("-" * 20 + "\n")
    
    # --- Step 4: Compare with Experimental Value ---
    ratio = tau_th / tau_exp
    
    print("--- Final Result ---")
    print(f"The ratio of theoretical to experimental lifetime is:")
    print(f"Ratio = τ_th / τ_exp = {tau_th:.3e} s / {tau_exp:.3e} s")
    print(f"Ratio ≈ {ratio:.3f}")
    
    # Determine which answer choice is closest
    choices = {'A': 30, 'B': 2, 'C': 0.1, 'D': 1, 'E': 100, 'F': 0.5, 'G': 10, 'H': 5}
    closest_choice = min(choices.keys(), key=lambda k: abs(choices[k] - ratio))
    print(f"\nThis multiple is closest to the value in answer choice '{closest_choice}': {choices[closest_choice]}.")

# Run the calculation
calculate_lifetime_and_ratio()