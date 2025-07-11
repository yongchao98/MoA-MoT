import math

def evaluate_quantum_correction_3d():
    """
    Calculates the quantum correction to conductivity (Weak Localization) for a
    bulk (3D) semiconductor and prints the detailed calculation.

    The correction is due to the constructive interference of time-reversed
    electron paths, enhancing backscattering and reducing conductivity. It is
    proportional to the probability of an electron returning to its origin.

    The final formula used is:
    Δσ = - (e² / (2π²ħ)) * (1/l_e - 1/L_φ)
    """

    # --- Physical Constants ---
    e = 1.602176634e-19  # Elementary charge in Coulombs (C)
    hbar = 1.054571817e-34 # Reduced Planck constant in Joule-seconds (J·s)
    pi = math.pi

    # --- Example Parameters for a Disordered Semiconductor ---
    # These values are illustrative. In a real scenario, they would be determined
    # from experimental measurements.
    # We must have L_phi > l_e for the correction to be physically meaningful.
    l_e = 10e-9   # Elastic mean free path in meters (m), e.g., 10 nm
    L_phi = 100e-9 # Phase coherence length in meters (m), e.g., 100 nm

    # --- Calculation ---
    # 1. Calculate the prefactor C = e² / (2π²ħ)
    prefactor = (e**2) / (2 * pi**2 * hbar)

    # 2. Calculate the length-dependent term (1/l_e - 1/L_phi)
    length_term = (1/l_e - 1/L_phi)

    # 3. Calculate the final conductivity correction Δσ
    delta_sigma = -prefactor * length_term

    # --- Output the detailed evaluation ---
    print("Evaluation of Quantum Correction to Conductivity (Weak Localization in 3D)")
    print("=======================================================================")
    print("The formula for the correction Δσ is derived from the integrated probability of electron return:")
    print("Δσ = - (e² / (2π²ħ)) * (1/l_e - 1/L_φ)\n")

    print("--- Used Values ---")
    print(f"e   (elementary charge)      = {e:.5e} C")
    print(f"ħ   (reduced Planck constant)  = {hbar:.5e} J·s")
    print(f"π   (pi)                     = {pi:.5f}")
    print(f"l_e (elastic mean free path) = {l_e:.1e} m")
    print(f"L_φ (phase coherence length) = {L_phi:.1e} m\n")

    print("--- Step-by-Step Calculation ---")
    print("1. Calculate the prefactor constant C = e² / (2π²ħ):")
    print(f"   C = ({e:.5e})² / (2 * {pi:.5f}² * {hbar:.5e})")
    print(f"   C = {prefactor:.5e} S (Siemens)\n")

    print("2. Calculate the length-dependent term ΔL⁻¹ = (1/l_e - 1/L_φ):")
    print(f"   ΔL⁻¹ = (1 / {l_e:.1e} - 1 / {L_phi:.1e})")
    print(f"   ΔL⁻¹ = {length_term:.3e} m⁻¹\n")

    print("3. Calculate the final correction Δσ = -C * ΔL⁻¹:")
    print(f"   Δσ = -({prefactor:.5e} S) * ({length_term:.3e} m⁻¹)")
    print(f"   Δσ = {delta_sigma:.2f} S/m\n")

    print("This negative value represents a decrease in the material's conductivity.")
    
    # Return the final numerical result for the 'answer' part
    return delta_sigma

# Execute the function to print the evaluation
final_result = evaluate_quantum_correction_3d()
# The final result is printed by the function, this captures the numeric value for the answer format.
# The format "<<<value>>>" should be at the very end of the entire response.
# final_answer_string = f"<<<{final_result:.2f}>>>"
# Let's check how this will be rendered. I should just output the string at the end.
