import math

def calculate_conductivity_correction():
    """
    Calculates the quantum correction to conductivity for an electron in a bulk semiconductor
    based on the theory of weak localization.
    """
    # Physical constants
    e = 1.602176634e-19  # Elementary charge in Coulombs
    hbar = 1.054571817e-34 # Reduced Planck constant in J*s
    pi = math.pi

    # Assumed typical parameters for a semiconductor (in meters)
    l_e = 10e-9  # Mean free path (10 nm)
    L_phi = 100e-9 # Phase coherence length (100 nm)

    # Calculate the prefactor C = e^2 / (2 * pi^2 * hbar)
    # The units of this prefactor are Siemens (S), which is a unit of conductance.
    prefactor = (e**2) / (2 * pi**2 * hbar)

    # Calculate the length-dependent term (1/l_e - 1/L_phi)
    # The units of this term are m^-1.
    inv_l_e = 1 / l_e
    inv_L_phi = 1 / L_phi
    length_term = inv_l_e - inv_L_phi

    # Calculate the final conductivity correction Δσ = -C * (1/l_e - 1/L_phi)
    # The units are S * m^-1, which is the correct unit for conductivity.
    delta_sigma = -prefactor * length_term

    # Print the evaluation step-by-step
    print("Evaluating the quantum correction to conductivity Δσ in 3D.")
    print("Formula: Δσ = - (e² / (2π²ħ)) * (1/lₑ - 1/Lᵩ)\n")
    print("Using the following values:")
    print(f"  e (elementary charge) = {e:.4e} C")
    print(f"  ħ (reduced Planck const) = {hbar:.4e} J·s")
    print(f"  lₑ (mean free path) = {l_e:.2e} m")
    print(f"  Lᵩ (phase coherence length) = {L_phi:.2e} m\n")

    print("Step 1: Calculate the prefactor (e² / (2π²ħ))")
    print(f"  Prefactor = ({e:.4e}² / (2 * π² * {hbar:.4e})) = {prefactor:.4e} S\n")
    
    print("Step 2: Calculate the length-dependent term (1/lₑ - 1/Lᵩ)")
    print(f"  Term = (1/{l_e:.2e} m - 1/{L_phi:.2e} m)")
    print(f"       = ({inv_l_e:.2e} m⁻¹ - {inv_L_phi:.2e} m⁻¹)")
    print(f"       = {length_term:.2e} m⁻¹\n")

    print("Step 3: Calculate the final conductivity correction Δσ")
    print(f"  Δσ = - ({prefactor:.4e} S) * ({length_term:.2e} m⁻¹)")
    print(f"  Δσ = {delta_sigma:.2f} S/m")

calculate_conductivity_correction()
<<< -1110.42 >>>