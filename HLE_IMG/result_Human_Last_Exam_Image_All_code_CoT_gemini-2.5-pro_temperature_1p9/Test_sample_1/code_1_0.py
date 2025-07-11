import math

def solve_reaction_kinetics():
    """
    Identifies the rate-determining step and calculates the rate constant for the given
    reaction based on the Energetic Span Model.
    """

    # Step 1: Define constants and energy values from the diagram.
    k_B = 1.380649e-23  # Boltzmann constant, J/K
    h = 6.62607015e-34   # Planck constant, J·s
    R = 8.31446          # Gas constant, J/(mol·K)
    T = 298.0            # Temperature, K

    # Define states with their Gibbs energy (G in kJ/mol) and sequence order
    intermediates = {
        'Cat':       {'G': 0,   'seq': 1},
        'Cat-GB':    {'G': -61, 'seq': 2},
        'R':         {'G': -55, 'seq': 3},
        'P5':        {'G': -22, 'seq': 5},
        'P':         {'G': -152,'seq': 7},
        'Cat-MPFA':  {'G': -142,'seq': 8}
    }
    transition_states = {
        'TS-add': {'G': 10, 'seq': 4},
        'TS-eli': {'G': 18, 'seq': 6}
    }

    # Step 2: Apply the Energetic Span Model to find TDI, TDTS, and the span.
    max_span = -float('inf')
    TDI_name = None
    TDTS_name = None

    for ts_name, ts_data in transition_states.items():
        for i_name, i_data in intermediates.items():
            # The model considers a span only if the intermediate precedes the transition state
            if i_data['seq'] < ts_data['seq']:
                span = ts_data['G'] - i_data['G']
                if span > max_span:
                    max_span = span
                    TDI_name = i_name
                    TDTS_name = ts_name

    delta_E_kJ = max_span
    TDI = intermediates[TDI_name]
    TDTS = transition_states[TDTS_name]

    # Map TDTS to its corresponding elementary step
    if TDTS_name == 'TS-eli':
        rds_description = "elimination of fluoride from P5 to form P"
    elif TDTS_name == 'TS-add':
        rds_description = "addition of water at state R to form P5"
    else:
        rds_description = "an unknown step"

    print("--- Identifying the Rate-Determining Step (RDS) ---")
    print("According to the Energetic Span Model, the rate of a catalytic cycle is determined by the largest energy difference between any transition state and any preceding intermediate.")
    print(f"\nThe analysis identifies:")
    print(f"1. The TOF-Determining Intermediate (TDI) is {TDI_name} with G = {TDI['G']} kJ/mol.")
    print(f"2. The TOF-Determining Transition State (TDTS) is {TDTS_name} with G = {TDTS['G']} kJ/mol.")
    print(f"\nThe rate-determining step is the one associated with the TDTS, which is the {rds_description}.")
    
    # Step 3: Calculate the Energetic Span (δE)
    delta_E_J = delta_E_kJ * 1000
    print("\n--- Calculating the Energetic Span (δE) ---")
    print(f"The energetic span (δE) is the Gibbs energy difference between the TDTS and the TDI:")
    print(f"δE = G({TDTS_name}) - G({TDI_name})")
    print(f"δE = {TDTS['G']} kJ/mol - ({TDI['G']}) kJ/mol = {delta_E_kJ} kJ/mol")
    print(f"δE = {delta_E_J:.0f} J/mol")

    # Step 4: Calculate the Rate Constant (k)
    print("\n--- Calculating the Rate Constant (k) ---")
    print("Using the Eyring equation: k = (k_B * T / h) * exp(-δE / (R * T))")
    
    pre_factor = (k_B * T) / h
    exponent = -delta_E_J / (R * T)
    k_s = pre_factor * math.exp(exponent)

    print("\nWith the values:")
    print(f"k_B * T / h = ({k_B:.4e} J/K * {T} K) / {h:.4e} J·s = {pre_factor:.3e} s⁻¹")
    print(f"-δE / (R * T) = -{delta_E_J:.0f} J/mol / ({R:.3f} J/(mol·K) * {T} K) = {exponent:.3f}")
    print(f"k = {pre_factor:.3e} s⁻¹ * exp({exponent:.3f}) = {k_s:.4f} s⁻¹")

    # Step 5: Convert Units and Format Output
    k_hr = k_s * 3600
    print(f"\nConverting rate constant units from s⁻¹ to h⁻¹:")
    print(f"k (h⁻¹) = k (s⁻¹) * 3600 s/h = {k_s:.4f} * 3600 = {k_hr:.2f} h⁻¹")

    # Round final answer to two significant figures
    final_answer = int(round(k_hr, 1 - int(math.floor(math.log10(abs(k_hr))))))
    
    print(f"\nThe calculated reaction rate constant to two significant figures is {final_answer} h⁻¹.")
    print(f"<<<{final_answer}>>>")

solve_reaction_kinetics()