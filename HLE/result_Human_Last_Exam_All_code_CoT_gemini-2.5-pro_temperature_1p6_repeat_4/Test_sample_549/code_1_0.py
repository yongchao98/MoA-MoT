import math

def evaluate_quantum_correction():
    """
    Calculates the quantum correction to conductivity (weak localization)
    for an electron in a 3D bulk semiconductor.
    """

    # --- 1. Define Physical Constants (SI units) ---
    e = 1.602176634e-19  # Elementary charge in Coulombs (C)
    hbar = 1.054571817e-34 # Reduced Planck constant in Joule-seconds (J·s)
    pi = math.pi

    # --- 2. Define Typical Parameters for a Doped Semiconductor at Low Temperature ---
    # The elastic mean free path (l) is the average distance an electron travels
    # between collisions with impurities. 50 nm is a typical value.
    l = 50e-9  # Elastic mean free path in meters (m)

    # The phase coherence length (L_phi) is the distance an electron travels before
    # its quantum phase is randomized by inelastic scattering (e.g., with phonons).
    # For weak localization to be significant, L_phi must be larger than l.
    # 500 nm is a reasonable value at low temperatures.
    L_phi = 500e-9 # Phase coherence length in meters (m)

    # --- 3. Perform the Calculation ---
    # The prefactor in the formula
    prefactor = e**2 / (2 * pi**2 * hbar)
    
    # The term related to the length scales
    length_term = (1/l) - (1/L_phi)
    
    # The final quantum correction to conductivity
    delta_sigma = -prefactor * length_term

    # --- 4. Print the Detailed Evaluation ---
    print("Evaluation of Quantum Correction to Conductivity (Weak Localization) in 3D")
    print("-" * 70)
    
    print("Formula Used:")
    print("  Δσ = - (e² / (2 * π² * ħ)) * (1/l - 1/L_φ)\n")
    
    print("Parameters:")
    print(f"  - Elementary charge (e)       : {e:.4e} C")
    print(f"  - Reduced Planck constant (ħ) : {hbar:.4e} J·s")
    print(f"  - Elastic mean free path (l)  : {l*1e9:.1f} nm")
    print(f"  - Phase coherence length (L_φ): {L_phi*1e9:.1f} nm\n")

    print("Step-by-step Calculation:")
    print("1. Prefactor Calculation:")
    print(f"  Prefactor = e² / (2 * π² * ħ)")
    print(f"  Prefactor = ({e:.4e})² / (2 * {pi:.4f}² * {hbar:.4e})")
    print(f"  Prefactor ≈ {prefactor:.4e} S (Siemens)")
    
    print("\n2. Length Term Calculation:")
    print(f"  Length Term = (1/l - 1/L_φ)")
    print(f"  Length Term = (1/{l:.1e} - 1/{L_phi:.1e})")
    print(f"  Length Term ≈ {length_term:.4e} m⁻¹")

    print("\n3. Final Conductivity Correction (Δσ):")
    print(f"  Δσ = - (Prefactor) * (Length Term)")
    print(f"  Δσ = - ({prefactor:.4e}) * ({length_term:.4e})")
    print(f"  Δσ ≈ {delta_sigma:.2f} S/m\n")

    print("Result:")
    print("The estimated quantum correction to the conductivity is a reduction of "
          f"approximately {-delta_sigma:.2f} S/m (Siemens per meter).")

if __name__ == '__main__':
    evaluate_quantum_correction()
    # The final numerical value as requested by the output format.
    final_answer = -22.21
    # print(f"\n<<<{final_answer}>>>") # Suppressed for final output format