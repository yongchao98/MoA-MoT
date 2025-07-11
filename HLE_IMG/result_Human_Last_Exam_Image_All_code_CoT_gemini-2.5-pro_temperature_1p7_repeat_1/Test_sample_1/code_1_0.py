import math

def calculate_rate_constant():
    """
    Identifies the rate-determining step and calculates the reaction rate constant
    based on the Energetic Span Model.
    """
    # Step 1: Identify the rate-determining step and energetic span (dG)
    # The Energetic Span Model states the rate is determined by the difference between the
    # highest transition state (TDTS) and the lowest intermediate (TDI) in the catalytic cycle.
    # The intermediates P (-152 kJ/mol) and Cat-MPFA (-142 kJ/mol) are deep thermodynamic sinks
    # representing product complexes. For the kinetic analysis of the turnover, we consider the
    # cycle preceding this irreversible sink.

    # Transition States (kJ/mol):
    G_TS_add = 10
    G_TS_eli = 18

    # Intermediates in the active cycle (kJ/mol):
    G_Cat = 0
    G_Cat_GB = -61
    G_R = -55
    G_P5 = -22

    # Identify the TDTS and TDI from the relevant species.
    # TDTS is the highest-energy TS: TS-eli
    G_TDTS = G_TS_eli  # +18 kJ/mol
    # TDI is the lowest-energy Intermediate in the active cycle: Cat-GB
    G_TDI = G_Cat_GB  # -61 kJ/mol
    
    # The rate-determining step (RDS) is the one containing the TDTS.
    # From the reaction diagram, TS-eli is the transition state for the P5 -> P step.
    print("1. Identifying the Rate-Determining Step (RDS) and Energetic Span:")
    print("   - Turnover-Determining Transition State (TDTS): TS-eli")
    print("   - Turnover-Determining Intermediate (TDI): Cat-GB")
    print("   - Rate-Determining Step (RDS): The elimination step P5 -> P.")

    # Calculate the energetic span in kJ/mol and J/mol
    dG_kj_mol = G_TDTS - G_TDI
    dG_j_mol = dG_kj_mol * 1000

    print(f"\n   Energetic Span (δE) = G(TDTS) - G(TDI) = {G_TDTS} kJ/mol - ({G_TDI} kJ/mol) = {dG_kj_mol} kJ/mol")

    # Step 2: Calculate the rate constant using the Eyring equation
    # k = (k_B * T / h) * exp(-dG / RT)

    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    h = 6.62607015e-34  # Planck constant (J*s)
    R = 8.31446         # Gas constant (J/(mol*K))
    
    # Reaction conditions
    T = 298  # Temperature (K)

    # Calculation
    pre_exponential_factor = (k_B * T) / h
    exponent = -dG_j_mol / (R * T)
    k_seconds = pre_exponential_factor * math.exp(exponent)
    
    # Convert from s^-1 to hours^-1
    k_hours = k_seconds * 3600
    
    # Round to two significant figures
    # For a number like 316.7, rounding to two significant figures gives 320.
    # We can achieve this by rounding to the nearest ten.
    k_hours_rounded = round(k_hours, -1)


    print("\n2. Calculating the Rate Constant (k):")
    print("   Using the Eyring equation: k = (k_B * T / h) * exp(-δE / (R * T))")
    print(f"   k = (({k_B:.4e} J/K * {T} K) / {h:.4e} J*s) * exp(-{dG_j_mol} J/mol / ({R:.4f} J/(mol*K) * {T} K))")
    print(f"\n   The calculated rate constant is {k_seconds:.2e} s^-1.")
    print(f"   In units of hours^-1, this is {k_hours:.1f} hours^-1.")
    print(f"\n   Rounded to two significant figures, the final reaction rate constant is {int(k_hours_rounded)} hours^-1.")
    
    return int(k_hours_rounded)

# Run the calculation and store the final answer
final_answer = calculate_rate_constant()