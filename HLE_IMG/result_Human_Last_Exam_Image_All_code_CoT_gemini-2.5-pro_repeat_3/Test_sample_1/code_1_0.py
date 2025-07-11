import math

def calculate_reaction_rate():
    """
    Identifies the rate-determining step using the Energetic Span Model
    and calculates the reaction rate constant.
    """
    # Step 1: Identify key species and energies from the provided diagram.
    # According to the Energetic Span Model, the rate of a catalytic cycle is determined by the
    # energy difference between the highest-energy transition state (TDTS) and the
    # lowest-energy intermediate (TDI).

    # From the Gibbs Energy profile:
    # The lowest-energy intermediate is P.
    G_TDI_kJ = -152  # Gibbs energy of intermediate P in kJ/mol
    # The highest-energy transition state is TS-eli.
    G_TDTS_kJ = 18   # Gibbs energy of transition state TS-eli in kJ/mol

    print("1. Identifying the Rate-Determining Step (RDS) and Energetic Span (ΔG‡)")
    print("--------------------------------------------------------------------------")
    print("According to the Energetic Span Model, the rate is determined by the highest energy transition state (TDTS) and the lowest energy intermediate (TDI).")
    print(f"From the diagram, the TDTS is TS-eli, with a Gibbs Energy of {G_TDTS_kJ} kJ/mol.")
    print(f"The TDI is the intermediate P, with a Gibbs Energy of {G_TDI_kJ} kJ/mol.")
    print("The rate-determining step is therefore the one involving the TDTS, which is the fluoride elimination step (P5 -> TS-eli -> P).")
    print("")

    # Step 2: Calculate the energetic span (ΔG‡).
    delta_G_ddagger_kJ = G_TDTS_kJ - G_TDI_kJ
    delta_G_ddagger_J = delta_G_ddagger_kJ * 1000  # Convert from kJ/mol to J/mol

    print("2. Calculating the Energetic Span (ΔG‡)")
    print("-----------------------------------------")
    print(f"The energetic span ΔG‡ is calculated as G(TDTS) - G(TDI).")
    print(f"ΔG‡ = {G_TDTS_kJ} kJ/mol - ({G_TDI_kJ} kJ/mol) = {delta_G_ddagger_kJ} kJ/mol")
    print(f"In Joules per mole, ΔG‡ = {delta_G_ddagger_J} J/mol.")
    print("")

    # Step 3: Calculate the rate constant k using the Eyring equation.
    # Define physical constants and temperature.
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34   # Planck constant in J*s
    R = 8.3144626        # Ideal gas constant in J/(mol*K)
    T = 298              # Temperature in K

    # Eyring equation: k = (k_B * T / h) * exp(-ΔG‡ / RT)
    pre_exponential_factor = (k_B * T) / h
    RT_term = R * T
    exponential_term = math.exp(-delta_G_ddagger_J / RT_term)
    k_seconds = pre_exponential_factor * exponential_term
    
    print("3. Calculating the Rate Constant (k)")
    print("---------------------------------------")
    print("The rate constant is calculated using the Eyring equation: k = (k_B * T / h) * exp(-ΔG‡ / RT)")
    print(f"Using the values:")
    print(f"k_B = {k_B:.6e} J/K")
    print(f"T = {T} K")
    print(f"h = {h:.6e} J*s")
    print(f"ΔG‡ = {delta_G_ddagger_J} J/mol")
    print(f"R = {R:.5f} J/(mol*K)")
    print("\nCalculation:")
    print(f"k = ({k_B:.3e} J/K * {T} K / {h:.3e} J*s) * exp(-{delta_G_ddagger_J} J/mol / ({R:.3f} J/(mol*K) * {T} K))")
    print(f"k = {k_seconds:.4e} s⁻¹")
    print("")

    # Step 4: Convert k to hours⁻¹ and format.
    k_hours = k_seconds * 3600

    print("4. Final Rate Constant in hours⁻¹")
    print("----------------------------------")
    print(f"To convert from s⁻¹ to hours⁻¹, multiply by 3600 s/hour.")
    print(f"k (hours⁻¹) = {k_seconds:.4e} s⁻¹ * 3600 s/hour = {k_hours:.4e} hours⁻¹")
    print("\n----------------------------------")
    print("Final Answer:")
    print(f"The reaction rate constant at 298 K, rounded to two significant figures, is {k_hours:.2g} hours⁻¹.")

if __name__ == '__main__':
    calculate_reaction_rate()