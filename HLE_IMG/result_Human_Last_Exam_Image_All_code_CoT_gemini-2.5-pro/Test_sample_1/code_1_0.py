import math

def solve_reaction_kinetics():
    """
    This script identifies the rate-determining step and calculates the rate constant
    for the given reaction based on the Energetic Span Model.
    """
    # Step 1: Define energies from the diagram (in kJ/mol)
    energies = {
        'Cat': 0, 
        'Cat-GB': -61, 
        'R': -55, 
        'TS-add': 10, 
        'P5': -22, 
        'TS-eli': 18, 
        'P': -152, 
        'Cat-MPFA': -142
    }

    # Step 2: Apply the Energetic Span Model
    # Identify the Turnover Determining Transition State (TDTS)
    ts_states = {k: v for k, v in energies.items() if 'TS' in k}
    tdts_name = max(ts_states, key=ts_states.get)
    g_tdts = ts_states[tdts_name]

    # Identify the Turnover Determining Intermediate (TDI)
    i_states = {k: v for k, v in energies.items() if 'TS' not in k}
    tdi_name = min(i_states, key=i_states.get)
    g_tdi = i_states[tdi_name]

    # The rate-determining step (RDS) is associated with the TDTS. 
    # From the diagram, TS-eli is the transition state for the P5 -> P step.
    rds_description = f"P5 -> P (via {tdts_name})"

    # Calculate the energetic span (ΔG‡) in kJ/mol
    delta_g_ddagger_kJ = g_tdts - g_tdi

    print("--- Energetic Span Model Analysis ---")
    print("\n1. Identification of the Rate-Determining Step (RDS):")
    print(f"   - The highest energy transition state (TDTS) is {tdts_name} at {g_tdts} kJ/mol.")
    print(f"   - The lowest energy intermediate (TDI) is {tdi_name} at {g_tdi} kJ/mol.")
    print(f"   - The rate-determining step is the one involving the TDTS, which is: {rds_description}")

    print("\n2. Calculation of the Energetic Span (Activation Energy, ΔG‡):")
    print(f"   - ΔG‡ = G({tdts_name}) - G({tdi_name})")
    print(f"   - ΔG‡ = {g_tdts} kJ/mol - ({g_tdi} kJ/mol) = {delta_g_ddagger_kJ:.0f} kJ/mol")

    # Step 3: Calculate the rate constant using the Eyring equation
    # Constants
    T = 298          # Temperature in Kelvin
    R = 8.314462     # Gas constant in J/(mol·K)
    k_B = 1.380649e-23 # Boltzmann constant in J/K
    h = 6.62607015e-34 # Planck constant in J·s
    
    # Convert ΔG‡ to J/mol
    delta_g_ddagger_J = delta_g_ddagger_kJ * 1000

    # Eyring equation
    k_seconds = (k_B * T / h) * math.exp(-delta_g_ddagger_J / (R * T))

    # Step 4: Convert rate constant to hours⁻¹
    k_hours = k_seconds * 3600

    print("\n3. Calculation of the Rate Constant (k):")
    print(f"   - Using the Eyring equation: k = (k_B * T / h) * exp(-ΔG‡ / (R * T))")
    print(f"   - The final equation with all values is:")
    print(f"     k = ({k_B:.6e} J/K * {T} K / {h:.6e} J·s) * exp(-{delta_g_ddagger_J:.0f} J/mol / ({R:.4f} J/(mol·K) * {T} K))")
    print(f"   - The calculated rate constant is {k_seconds:.2e} s⁻¹.")
    print(f"   - In hours⁻¹, k = {k_hours:.2e} hours⁻¹.")

    # Format final answer to two significant figures
    final_answer = f"{k_hours:.1e}"
    
    print("\n--- Final Answer ---")
    print(f"The rate-determining step is {rds_description}.")
    print(f"The reaction rate constant is {final_answer} hours⁻¹.")
    
    # Submission format
    print(f"<<<{final_answer}>>>")

solve_reaction_kinetics()