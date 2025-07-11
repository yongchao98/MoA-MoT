import math

def calculate_reaction_rate():
    """
    Identifies the RDS and calculates the rate constant for the given reaction
    based on the Energetic Span Model.
    """
    # Step 1: Identify TDTS, TDI, and RDS from the energy diagram.
    # The highest energy transition state (TDTS) is TS-eli.
    G_TDTS_kJ_mol = 18.0
    # The lowest energy intermediate (TDI) is P.
    G_TDI_kJ_mol = -152.0

    print("Identification of the Rate-Determining Step (RDS) using the Energetic Span Model:")
    print("The highest-energy transition state (TDTS) is TS-eli at +18 kJ/mol.")
    print("The lowest-energy intermediate (TDI) is P at -152 kJ/mol.")
    print("The rate-determining step is the one involving the TDTS, which is the elimination step (P5 -> P).")
    print("-" * 40)

    # Step 2: Calculate the energetic span (δE).
    delta_E_kJ_mol = G_TDTS_kJ_mol - G_TDI_kJ_mol

    # Step 3: Calculate the rate constant using the Eyring equation.
    # Define physical constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34   # Planck constant in J*s
    R = 8.3144626       # Ideal gas constant in J/(mol*K)
    T = 298.0           # Temperature in K

    # Convert energetic span from kJ/mol to J/mol
    delta_E_J_mol = delta_E_kJ_mol * 1000

    # Calculate the pre-exponential and exponential terms
    pre_exponential_factor = (k_B * T) / h
    exponent = -delta_E_J_mol / (R * T)
    k_seconds = pre_exponential_factor * math.exp(exponent)

    # Step 4: Convert the rate constant to hours^-1
    k_hours = k_seconds * 3600

    # Step 5: Print the detailed calculation and the final result.
    print("Calculation of the rate constant (k):")
    print("\nEquation: k = (k_B * T / h) * exp(-δE / (R * T))")
    
    print("\nValues used:")
    print(f"Energetic Span (δE) = G(TDTS) - G(TDI)")
    print(f"δE = {G_TDTS_kJ_mol} kJ/mol - ({G_TDI_kJ_mol} kJ/mol) = {delta_E_kJ_mol} kJ/mol = {delta_E_J_mol} J/mol")
    print(f"T = {T} K")
    print(f"k_B = {k_B:.6e} J/K")
    print(f"h = {h:.6e} J*s")
    print(f"R = {R:.5f} J/(mol*K)")
    
    print("\nFinal Equation with numbers:")
    print(f"k = (({k_B:.6e} * {T}) / {h:.6e}) * exp(-{delta_E_J_mol} / ({R:.5f} * {T}))")

    print("\nResult:")
    print(f"The calculated rate constant is {k_hours:.2g} hours^-1.")


calculate_reaction_rate()
<<<2.9e-14>>>