import math

def evaluate_conductivity_correction():
    """
    Evaluates the quantum correction to conductivity in a bulk semiconductor.
    """
    # 1. Explanation of the physical model and formula
    print("This script evaluates the quantum correction to the conductivity (Δσ) for an electron in a bulk (3D) semiconductor.")
    print("The effect, known as weak localization, results from the constructive quantum interference of an electron traveling along a closed-loop path with its time-reversed counterpart.")
    print("This interference increases the electron's return probability, leading to an increase in resistance (a decrease in conductivity).\n")
    
    print("The derived formula for this correction is:")
    print("Δσ = - (e² / (π² * ħ)) * [1/l - (1/L_φ) * arctan(L_φ/l)]\n")
    print("Where:")
    print("  Δσ  : The change in conductivity (in Siemens/meter)")
    print("  e   : The elementary charge (1.602 x 10⁻¹⁹ C)")
    print("  ħ   : The reduced Planck constant (1.054 x 10⁻³⁴ J·s)")
    print("  π   : The mathematical constant Pi")
    print("  l   : The elastic mean free path of the electron (in meters)")
    print("  L_φ : The phase coherence length (in meters)\n")
    print("The correction Δσ is negative, confirming that weak localization reduces the conductivity.")
    print("-" * 50)

    # 2. Sample Calculation
    print("Sample Calculation for a Doped Semiconductor at Low Temperature\n")

    # Physical constants
    e = 1.60217663e-19      # Coulombs
    hbar = 1.05457182e-34   # Joule-seconds
    pi = math.pi

    # Typical experimental parameters for a semiconductor
    l = 100e-9      # Elastic mean free path (100 nm)
    L_phi = 1000e-9 # Phase coherence length (1 µm or 1000 nm)

    print("Using the following parameters:")
    print(f"  e       = {e:.4e} C")
    print(f"  ħ       = {hbar:.4e} J·s")
    print(f"  π       = {pi:.4f}")
    print(f"  l       = {l * 1e9:.0f} nm ({l:.1e} m)")
    print(f"  L_φ     = {L_phi * 1e9:.0f} nm ({L_phi:.1e} m)\n")

    # Perform the calculation step-by-step
    # Step 1: Calculate the prefactor constant
    prefactor = e**2 / (pi**2 * hbar)
    
    # Step 2: Calculate the term in the brackets, which depends on the lengths
    length_term = (1 / l) - (1 / L_phi) * math.atan(L_phi / l)
    
    # Step 3: Calculate the final conductivity correction
    delta_sigma = -prefactor * length_term

    # Print the detailed calculation
    print("Calculation Steps:")
    
    # Print numbers in the final equation
    print("1. Prefactor constant, (e² / (π² * ħ)):")
    print(f"   ({e:.4e}² / ({pi:.4f}² * {hbar:.4e})) = {prefactor:.4e} S (Siemens)")

    print("\n2. Length-dependent term, [1/l - (1/L_φ) * arctan(L_φ/l)]:")
    print(f"   [1/{l:.1e} - (1/{L_phi:.1e}) * arctan({L_phi/l:.1f})] = {length_term:.4e} m⁻¹")
    
    print("\n3. Final Quantum Correction, Δσ = - (Prefactor) * (Length Term):")
    print(f"   Δσ = -({prefactor:.4e}) * ({length_term:.4e})")
    
    print("-" * 50)
    print(f"The evaluated quantum correction to the conductivity is: Δσ = {delta_sigma:.2f} S/m")
    print("-" * 50)

# Execute the function
evaluate_conductivity_correction()