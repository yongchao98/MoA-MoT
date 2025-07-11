import math

def calculate_conductivity_correction():
    """
    Evaluates the quantum correction to conductivity (weak localization)
    for an electron in a 3D semiconductor.
    """
    # Physical constants
    e = 1.60217663e-19  # Elementary charge in Coulombs
    hbar = 1.05457182e-34 # Reduced Planck constant in J*s
    pi = math.pi

    # Assumed parameters for a typical bulk semiconductor at low temperature
    # The elastic mean free path (l_e)
    l_e = 100e-9  # 100 nm in meters
    # The phase coherence length (L_phi)
    L_phi = 1000e-9 # 1000 nm (1 um) in meters

    # The formula for the 3D quantum correction to conductivity is:
    # Δσ = - (e^2 / (2 * π^2 * ħ)) * (1/l_e - 1/L_phi)

    # Step 1: Calculate the prefactor C = e^2 / (2 * π^2 * ħ)
    # The units of this prefactor are Siemens (S).
    prefactor = e**2 / (2 * pi**2 * hbar)

    # Step 2: Calculate the length-dependent term (1/l_e - 1/L_phi)
    # The units of this term are m^-1.
    length_term = (1 / l_e) - (1 / L_phi)

    # Step 3: Calculate the final quantum correction Δσ
    # The units are S * m^-1 = S/m.
    delta_sigma = -prefactor * length_term
    
    print("Evaluating the quantum correction to conductivity (Weak Localization) in 3D.")
    print("Formula: Δσ = - (e^2 / (2 * π^2 * ħ)) * (1/l_e - 1/L_phi)\n")
    print("Assumed parameters:")
    print(f"  - Elastic mean free path (l_e): {l_e:.0e} m")
    print(f"  - Phase coherence length (L_phi): {L_phi:.0e} m\n")

    print("Calculation steps:")
    # Print the equation with numerical values
    # First, the prefactor part
    print(f"1. Prefactor C = e^2 / (2 * π^2 * ħ):")
    print(f"   C = ({e:.4e} C)^2 / (2 * {pi:.4f}^2 * {hbar:.4e} J*s) = {prefactor:.4e} S")

    # Second, the length-dependent part
    print(f"\n2. Length-dependent term = (1/l_e - 1/L_phi):")
    print(f"   Term = (1/{l_e:.1e} m - 1/{L_phi:.1e} m) = {length_term:.4e} m^-1")

    # Finally, the full calculation
    print("\n3. Final result Δσ = -C * (Term):")
    print(f"   Δσ = -({prefactor:.4e} S) * ({length_term:.4e} m^-1)")
    print(f"   Δσ = {delta_sigma:.2f} S/m\n")
    
    # Return the final numeric value for the answer tag
    return delta_sigma

# Execute the calculation and print the final answer in the required format
final_answer = calculate_conductivity_correction()
print(f"<<<{final_answer:.2f}>>>")
