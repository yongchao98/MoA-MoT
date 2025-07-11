import math

def calculate_reaction_rate():
    """
    This function identifies the rate-determining step and calculates the 
    reaction rate constant based on the provided energy profile.
    """

    # 1. Define constants and Gibbs free energies from the diagram
    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck constant in J*s
    R = 8.3144626     # Ideal gas constant in J/(mol*K)
    T = 298             # Temperature in Kelvin

    # Energies from the diagram in kJ/mol
    G_R = -55
    G_TS_add = 10
    G_P5 = -22
    G_TS_eli = 18

    # 2. Identify the Rate-Determining Step (RDS) by comparing activation barriers
    print("--- Step 1: Identifying the Rate-Determining Step ---")
    
    # Calculate activation barrier for the addition step (R -> P5)
    delta_G_add = G_TS_add - G_R
    print(f"Activation energy for addition (R -> P5) = G(TS-add) - G(R) = {G_TS_add} - ({G_R}) = {delta_G_add} kJ/mol")

    # Calculate activation barrier for the elimination step (P5 -> P)
    delta_G_eli = G_TS_eli - G_P5
    print(f"Activation energy for elimination (P5 -> P) = G(TS-eli) - G(P5) = {G_TS_eli} - ({G_P5}) = {delta_G_eli} kJ/mol")

    # Determine the RDS
    if delta_G_add > delta_G_eli:
        rds_step = "Addition (R -> P5)"
        delta_G_activation_kJ = delta_G_add
    else:
        rds_step = "Elimination (P5 -> P)"
        delta_G_activation_kJ = delta_G_eli

    print(f"\nComparing the barriers, the highest is {delta_G_activation_kJ} kJ/mol.")
    print(f"The Rate-Determining Step (RDS) is: {rds_step}")

    # 3. Calculate the rate constant using the Eyring equation
    print("\n--- Step 2: Calculating the Reaction Rate Constant (k) ---")
    
    # Convert activation energy to J/mol
    delta_G_activation_J = delta_G_activation_kJ * 1000
    
    print(f"The activation energy for the rate-determining step is ΔG‡ = {delta_G_activation_kJ} kJ/mol.")
    print(f"Using the Eyring equation: k = (k_B * T / h) * exp(-ΔG‡ / (R * T))")
    
    # Eyring equation calculation
    pre_exponential_factor = (k_B * T) / h
    exponent = -delta_G_activation_J / (R * T)
    k_s = pre_exponential_factor * math.exp(exponent)
    
    print(f"k = (({k_B:.4e} J/K * {T} K) / {h:.4e} J*s) * exp(-{delta_G_activation_J} J/mol / ({R:.3f} J/(mol*K) * {T} K))")
    print(f"k = ({pre_exponential_factor:.4e} s^-1) * exp({exponent:.4f})")
    print(f"k = {k_s:.4f} s^-1")
    
    # 4. Convert units and format the final answer
    k_hr = k_s * 3600
    
    print("\n--- Step 3: Final Answer ---")
    print(f"The rate constant in hours^-1 is: {k_s:.4f} s^-1 * 3600 s/hr = {k_hr:.4f} hr^-1")
    print(f"Formatted to two significant figures, the rate constant is: {k_hr:.2g} hours^-1")

if __name__ == '__main__':
    calculate_reaction_rate()
<<<8.9e+04>>>