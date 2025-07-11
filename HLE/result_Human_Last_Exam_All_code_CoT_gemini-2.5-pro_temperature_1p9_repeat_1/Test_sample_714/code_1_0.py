import numpy as np

def calculate_lifetime():
    """
    Calculates the theoretical lifetime of the Sodium 3p state using a hydrogenic model,
    and compares it to the experimental value.
    """

    # --- Step 1: Define Constants and Parameters ---
    # Physical constants
    alpha = 7.297e-3  # Fine-structure constant
    c = 3.00e8       # Speed of light in m/s
    a0 = 5.29e-11    # Bohr radius in m
    
    # Given parameters
    lambda_nm = 589  # Wavelength in nm
    lambda_m = lambda_nm * 1e-9 # Wavelength in m
    tau_exp_ns = 16.2 # Experimental lifetime in ns
    tau_exp_s = tau_exp_ns * 1e-9 # Experimental lifetime in s
    
    # Effective nuclear charge for the valence electron in Sodium.
    # In the hydrogenic model approximation for alkali metals, Z_eff=1 is used.
    Z_eff = 1.0

    print("--- Plan ---")
    print("1. Calculate the transition's angular frequency (omega) from its wavelength.")
    print("2. Calculate the theoretical transition rate (A) using the formula for hydrogenic atoms.")
    print("   Note: A crucial assumption is using an effective nuclear charge Z_eff=1 for the valence electron.")
    print("3. Calculate the theoretical lifetime (tau) as 1/A.")
    print("4. Compare the theoretical lifetime to the experimental lifetime.")
    print("-" * 20)

    # --- Step 2: Calculate Angular Frequency (omega) ---
    omega = 2 * np.pi * c / lambda_m
    print(f"--- Calculating Angular Frequency (ω) ---")
    print(f"ω = 2 * π * c / λ")
    print(f"ω = 2 * {np.pi:.4f} * {c:.2e} m/s / {lambda_m:.2e} m")
    print(f"ω = {omega:.4e} rad/s\n")

    # --- Step 3: Calculate the Transition Rate (A) ---
    # For a 3p -> 3s transition, the squared radial integral |<3s|r|3p>|^2 = 162 * (a0/Z_eff)^2
    # The formula for the rate A simplifies to: A = 72 * alpha * omega^3 * a0^2 / (c^2 * Z_eff^2)
    print(f"--- Calculating Theoretical Transition Rate (A) ---")
    print(f"The formula for the rate of a 3p->3s transition is:")
    print(f"A = (72 * α * ω³ * a₀²) / (c² * Z_eff²)\n")

    print(f"Plugging in the values (with Z_eff = {Z_eff}):")
    term1_str = "72"
    term2_str = f"{alpha:.4e}"
    term3_str = f"({omega:.4e} rad/s)³"
    term4_str = f"({a0:.2e} m)²"
    denom1_str = f"({c:.2e} m/s)²"
    denom2_str = f"({Z_eff})²"
    print(f"A = ({term1_str} * {term2_str} * {term3_str} * {term4_str}) / ({denom1_str} * {denom2_str})")
    
    A = (72 * alpha * (omega**3) * (a0**2)) / (c**2 * Z_eff**2)
    
    numerator_val = 72 * alpha * (omega**3) * (a0**2)
    denominator_val = (c**2 * Z_eff**2)
    
    print(f"A = {numerator_val:.4e} / {denominator_val:.4e}")
    print(f"A = {A:.4e} s⁻¹\n")

    # --- Step 4: Calculate Theoretical Lifetime (tau_th) ---
    tau_th_s = 1 / A
    tau_th_ns = tau_th_s * 1e9
    print(f"--- Calculating Theoretical Lifetime (τ_th) ---")
    print(f"τ_th = 1 / A = 1 / {A:.4e} s⁻¹ = {tau_th_s:.4e} s")
    print(f"τ_th = {tau_th_ns:.2f} ns\n")
    
    # --- Step 5: Compare Theoretical and Experimental Lifetimes ---
    ratio = tau_th_s / tau_exp_s
    print(f"--- Comparison with Experiment ---")
    print(f"Experimental lifetime τ_exp = {tau_exp_ns} ns")
    print(f"Theoretical lifetime τ_th = {tau_th_ns:.2f} ns")
    print(f"The ratio is τ_th / τ_exp = {tau_th_ns:.2f} / {tau_exp_ns} = {ratio:.2f}\n")
    print(f"The calculated theoretical lifetime is approximately {ratio:.1f} times the experimental lifetime.")
    print("This corresponds to 'one tenth as long as experiment'.")

calculate_lifetime()
<<<C>>>