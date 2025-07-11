import math

def solve_reaction_kinetics():
    """
    Identifies the rate-determining step using the Energetic Span Model
    and calculates the reaction rate constant for the given reaction.
    """
    # Step 1: Define energies from the reaction profile (in kJ/mol)
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

    # Define the sequence of species in the reaction path
    # This helps determine which intermediates precede a transition state.
    reaction_path = ['Cat', 'Cat-GB', 'R', 'TS-add', 'P5', 'TS-eli', 'P', 'Cat-MPFA']
    
    # Identify intermediates and transition states
    intermediates = {k: v for k, v in energies.items() if 'TS' not in k}
    transition_states = {k: v for k, v in energies.items() if 'TS' in k}
    
    print("--- Step 1: Applying the Energetic Span Model ---")
    
    max_span = -1
    tdi = None
    tdts = None
    
    # Iterate through each transition state to find the energetic span
    for ts_name, ts_energy in transition_states.items():
        ts_index = reaction_path.index(ts_name)
        preceding_intermediates = {k: v for k, v in intermediates.items() if reaction_path.index(k) < ts_index}
        
        # Find the lowest-energy intermediate preceding this TS
        if preceding_intermediates:
            min_energy_i_name = min(preceding_intermediates, key=preceding_intermediates.get)
            min_energy_i_value = preceding_intermediates[min_energy_i_name]
            
            # Calculate the span
            current_span = ts_energy - min_energy_i_value
            print(f"Calculating span for {ts_name}:")
            print(f"  Preceding intermediates: {list(preceding_intermediates.keys())}")
            print(f"  Lowest preceding intermediate (TDI candidate): {min_energy_i_name} (G = {min_energy_i_value} kJ/mol)")
            print(f"  Transition state (TDTS candidate): {ts_name} (G = {ts_energy} kJ/mol)")
            print(f"  Calculated Span = {ts_energy} - ({min_energy_i_value}) = {current_span} kJ/mol\n")

            if current_span > max_span:
                max_span = current_span
                tdi = (min_energy_i_name, min_energy_i_value)
                tdts = (ts_name, ts_energy)

    energetic_span_dE = max_span # in kJ/mol
    
    # Determine the RDS from the TDTS
    # TS-eli is the transition state between P5 and P.
    rds = "P5 -> P"
    
    print("--- Energetic Span Model Results ---")
    print(f"The TOF-Determining Intermediate (TDI) is {tdi[0]} with a Gibbs energy of {tdi[1]} kJ/mol.")
    print(f"The TOF-Determining Transition State (TDTS) is {tdts[0]} with a Gibbs energy of {tdts[1]} kJ/mol.")
    print(f"The Energetic Span (δE) is G(TDTS) - G(TDI) = {tdts[1]} - ({tdi[1]}) = {energetic_span_dE} kJ/mol.")
    print(f"Therefore, the rate-determining step is the elimination step: {rds}\n")
    
    # --- Step 2: Calculate the rate constant k using the Eyring equation ---
    print("--- Step 2: Calculating the Rate Constant (k) ---")
    
    # Constants
    k_B = 1.380649e-23 # Boltzmann constant, J/K
    h = 6.62607015e-34 # Planck's constant, J*s
    R = 8.314462618   # Ideal gas constant, J/(mol*K)
    T = 298            # Temperature, K
    
    # Convert energetic span from kJ/mol to J/mol
    dE_J = energetic_span_dE * 1000
    
    print("Using the Eyring equation: k = (k_B * T / h) * exp(-δE / (R * T))")
    print("Plugging in the values:")
    print(f"k_B = {k_B:.3e} J/K")
    print(f"T = {T} K")
    print(f"h = {h:.3e} J*s")
    print(f"δE = {dE_J} J/mol")
    print(f"R = {R:.4f} J/(mol*K)\n")

    # Calculation
    pre_exponential_factor = (k_B * T) / h
    exponent = -dE_J / (R * T)
    rate_constant_s = pre_exponential_factor * math.exp(exponent)
    
    print("Calculation steps:")
    print(f"k = ({k_B:.3e} * {T} / {h:.3e}) * exp(-{dE_J} / ({R:.4f} * {T}))")
    print(f"k = ({pre_exponential_factor:.3e}) * exp({exponent:.3f})")
    print(f"k = {rate_constant_s:.4f} s⁻¹\n")

    # Convert to hours^-1
    rate_constant_hr = rate_constant_s * 3600
    
    print("Converting the rate constant to hours⁻¹:")
    print(f"k (hours⁻¹) = k (s⁻¹) * 3600 s/hr")
    print(f"k (hours⁻¹) = {rate_constant_s:.4f} * 3600 = {rate_constant_hr:.2f} hours⁻¹\n")
    
    # Round to two significant figures
    rate_constant_hr_2sf = float(f"{rate_constant_hr:.1e}")
    
    print("--- Final Answer ---")
    print(f"The calculated reaction rate constant is {rate_constant_hr_2sf:.1e} hours⁻¹ (rounded to two significant figures).")
    
    # Final answer in specified format
    print(f"<<<{rate_constant_hr_2sf:.1e}>>>")

solve_reaction_kinetics()