import math

def calculate_rate_constant():
    """
    Identifies the rate-determining step based on the Energetic Span Model
    and calculates the reaction rate constant.
    """
    # --- Step 1: Define Constants and Input Energies ---
    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34 # Planck constant in J·s
    R = 8.3144626      # Ideal gas constant in J/(mol·K)
    T = 298            # Temperature in K

    # Energies from the provided diagram (in kJ/mol)
    # The TDTS is the highest-energy transition state in the cycle.
    G_TDTS = 18  # Energy of TS-eli
    # The TDI is the lowest-energy intermediate in the cycle.
    G_TDI = -152 # Energy of P

    # --- Step 2: Identify RDS and Calculate Energetic Span (delta_E) ---
    print("--- Rate-Determining Step Identification (Energetic Span Model) ---")
    print("The model defines the rate-determining span by the highest energy transition state (TDTS) and the lowest energy intermediate (TDI) in the cycle.")
    print(f"The TDTS is TS-eli with a Gibbs Free Energy of {G_TDTS} kJ/mol.")
    print(f"The TDI is the intermediate P with a Gibbs Free Energy of {G_TDI} kJ/mol.\n")

    delta_E_kJ = G_TDTS - G_TDI
    print("--- Calculation of the Energetic Span (delta_E) ---")
    print("The formula for the energetic span is: delta_E = G(TDTS) - G(TDI)")
    print(f"delta_E = {G_TDTS} kJ/mol - ({G_TDI} kJ/mol) = {delta_E_kJ} kJ/mol\n")

    # Convert delta_E to J/mol for the Eyring equation
    delta_E_J = delta_E_kJ * 1000

    # --- Step 3: Calculate the Rate Constant (k) ---
    # The rate constant 'k' is calculated using the Eyring equation with the energetic span.
    # Formula: k = (k_B * T / h) * exp(-delta_E / (R * T))

    print("--- Calculation of the Rate Constant (k) ---")
    print("The formula for the rate constant is: k = (k_B * T / h) * exp(-delta_E / (R * T))")
    # Printing the equation with values plugged in
    print("Plugging in the values:")
    print(f"k = (({k_B:.3e} J/K * {T} K) / {h:.3e} J·s) * exp(-{delta_E_J:.0f} J/mol / ({R:.3f} J/(mol·K) * {T} K))\n")

    # Perform the calculation
    k_s = (k_B * T / h) * math.exp(-delta_E_J / (R * T))

    # Convert rate from s^-1 to h^-1 (1 hour = 3600 seconds)
    k_h = k_s * 3600

    # --- Step 4: Final Output ---
    print("--- Final Result ---")
    print(f"The reaction rate constant is {k_h:.1e} hours^-1.")


if __name__ == '__main__':
    calculate_rate_constant()