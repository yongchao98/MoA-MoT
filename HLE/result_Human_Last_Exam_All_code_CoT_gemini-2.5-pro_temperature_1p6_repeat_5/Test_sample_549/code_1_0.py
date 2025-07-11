import math

def evaluate_conductivity_correction():
    """
    Calculates the quantum correction to conductivity for an electron
    in a bulk (3D) semiconductor due to weak localization.
    """
    # --- 1. Physical Constants (in SI units) ---
    e = 1.602176634e-19  # Elementary charge in Coulombs (C)
    hbar = 1.054571817e-34 # Reduced Planck constant in Joule-seconds (J·s)
    pi = math.pi

    # --- 2. Input Parameters (Typical values for a doped semiconductor) ---
    # These values can be changed to model different materials or conditions.
    fermi_velocity_vf = 3.0e5       # Fermi velocity in meters/second (m/s)
    elastic_time_tau_e = 1.0e-13    # Elastic scattering time in seconds (s)
    phase_coherence_time_tau_phi = 1.0e-11 # Phase coherence time in seconds (s)

    # --- 3. Calculate Intermediate Physical Quantities ---
    # Diffusion constant in 3D
    diffusion_constant_D = (fermi_velocity_vf**2 * elastic_time_tau_e) / 3

    # Elastic mean free path
    l_e = fermi_velocity_vf * elastic_time_tau_e

    # Phase coherence length
    L_phi = math.sqrt(diffusion_constant_D * phase_coherence_time_tau_phi)

    # --- 4. Calculate the Quantum Correction to Conductivity (Δσ) ---
    # Ensure L_phi > l_e, otherwise the model is not applicable
    if L_phi <= l_e:
        print("Warning: Phase coherence length is not greater than the mean free path.")
        print("The weak localization model may not be applicable.")
        delta_sigma = 0
    else:
        # Pre-factor in the formula
        prefactor = e**2 / (2 * pi**2 * hbar)
        # The length-dependent term
        length_term = (1/l_e) - (1/L_phi)
        # The final correction is negative, indicating a reduction in conductivity
        delta_sigma = -prefactor * length_term


    # --- 5. Print the Results ---
    print("--- Evaluation of Quantum Correction to Conductivity (Weak Localization in 3D) ---")
    print("\n[Input Parameters]")
    print(f"Fermi Velocity (v_F): {fermi_velocity_vf:.2e} m/s")
    print(f"Elastic Scattering Time (τ_e): {elastic_time_tau_e:.2e} s")
    print(f"Phase Coherence Time (τ_φ): {phase_coherence_time_tau_phi:.2e} s")

    print("\n[Intermediate Calculations]")
    print(f"Diffusion Constant (D): {diffusion_constant_D:.3e} m²/s")
    print(f"Elastic Mean Free Path (l_e): {l_e:.3e} m")
    print(f"Phase Coherence Length (L_φ): {L_phi:.3e} m")

    print("\n[Final Equation and Calculation]")
    print(f"Formula: Δσ = - (e² / (2 * π² * ħ)) * (1/lₑ - 1/Lφ)")
    # Using single-letter variables for clarity in the equation string
    equation_str = (
        f"Δσ = - (({e:.3e} C)² / (2 * ({pi:.4f})² * {hbar:.3e} J·s)) * (1 / {l_e:.3e} m - 1 / {L_phi:.3e} m)\n"
        f"Δσ = - ({prefactor:.3e} C²/(J·s)) * ({length_term:.3e} m⁻¹)"
    )
    print(equation_str)

    print("\n[Result]")
    print(f"The quantum correction to conductivity (Δσ) is: {delta_sigma:.4f} S/m")

    # Return the final value for the answer tag
    return delta_sigma

# Execute the function and capture the result
final_answer = evaluate_conductivity_correction()
print(f"\n<<<_EOF_{final_answer}>>>")