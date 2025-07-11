import math

def solve_reaction_kinetics():
    """
    This function identifies the rate-determining step, calculates the energetic span,
    and then computes the reaction rate constant based on the provided energy profile.
    """

    # --- Step 1: Define constants and given data ---
    # Physical constants
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34   # Planck constant in J*s
    R = 8.3144626      # Ideal gas constant in J/(mol*K)

    # Reaction conditions from the problem
    T = 298  # Temperature in Kelvin

    # Gibbs Free Energies from the diagram (in kJ/mol)
    energies = {
        'Cat': 0, 'Cat-GB': -61, 'R': -55, 'TS-add': 10,
        'P5': -22, 'TS-eli': 18, 'P': -152, 'Cat-MPFA': -142
    }
    
    # Define the sequence of species along the reaction coordinate
    reaction_path = ['Cat', 'Cat-GB', 'R', 'TS-add', 'P5', 'TS-eli', 'P', 'Cat-MPFA']

    # --- Step 2: Identify Rate-Determining Step and Energetic Span (δE) ---
    print("Step 1: Identify the rate-determining step and energetic span (δE) using the Energetic Span Model.")

    # Separate species into intermediates and transition states
    transition_states = {name: energy for name, energy in energies.items() if 'TS' in name}
    intermediates = {name: energy for name, energy in energies.items() if 'TS' not in name}

    # Find the highest-energy transition state (TDTS)
    tdts_name = max(transition_states, key=transition_states.get)
    tdts_energy = transition_states[tdts_name]

    # To find the rate-determining intermediate (TDI), we look for the most stable intermediate
    # that appears *before* the TDTS in the reaction pathway.
    tdts_index = reaction_path.index(tdts_name)
    preceding_species_names = reaction_path[:tdts_index]
    preceding_intermediates = {name: energies[name] for name in preceding_species_names if name in intermediates}
    
    # The TDI is the lowest energy state among these preceding intermediates.
    tdi_name = min(preceding_intermediates, key=preceding_intermediates.get)
    tdi_energy = preceding_intermediates[tdi_name]
    
    # The rate-determining step is the one associated with the TDTS.
    # Looking at the diagram, TS-eli is the transition state for the P5 -> P elimination step.
    print(f"The rate-determining step is the elimination step, which proceeds through the highest-energy transition state, {tdts_name}.")
    
    # Calculate the energetic span (δE)
    delta_E_kJ = tdts_energy - tdi_energy
    
    print("\n--- Energetic Span Calculation ---")
    print(f"Highest-energy Transition State (TDTS): {tdts_name} at G = {tdts_energy} kJ/mol")
    print(f"Rate-determining Intermediate (TDI): {tdi_name} at G = {tdi_energy} kJ/mol")
    print(f"The energetic span δE = G(TDTS) - G(TDI)")
    print(f"δE = {tdts_energy} kJ/mol - ({tdi_energy} kJ/mol) = {delta_E_kJ} kJ/mol")
    

    # --- Step 3: Calculate the reaction rate constant (k) ---
    print("\n\nStep 2: Calculate the reaction rate constant (k) using the Eyring equation.")

    # Convert energetic span from kJ/mol to J/mol for calculation
    delta_E_J = delta_E_kJ * 1000

    # Eyring equation: k = (k_B * T / h) * exp(-δE / (R * T))
    # Calculate k in units of s^-1
    k_s = (k_B * T / h) * math.exp(-delta_E_J / (R * T))

    # Convert rate constant k from s^-1 to hr^-1
    k_hr = k_s * 3600

    print("\n--- Rate Constant Calculation ---")
    print("The Eyring equation is: k = (k_B * T / h) * exp(-δE / (R * T))")
    print("Plugging in the values:")
    print(f"k = ({k_B:.5e} J/K * {T} K / {h:.5e} J*s) * exp(-{delta_E_J} J/mol / ({R:.5f} J/(mol*K) * {T} K))")
    
    # Calculate intermediate parts of the equation
    prefactor = k_B * T / h
    exponent_val = -delta_E_J / (R * T)
    print(f"k = ({prefactor:.3e} s⁻¹) * exp({exponent_val:.3f})")
    print(f"k = {k_s:.3f} s⁻¹")
    print(f"Converting to hr⁻¹: k = {k_s:.3f} s⁻¹ * 3600 s/hr = {k_hr:.2f} hr⁻¹")
    
    # Round the final answer to two significant figures
    # For a number like 319.4, rounding to 2 s.f. gives 320.
    k_hr_rounded = round(k_hr, -1)

    print(f"\nThe final calculated reaction rate constant, rounded to two significant figures, is {k_hr_rounded} hours⁻¹.")
    
    return k_hr_rounded

if __name__ == '__main__':
    final_answer = solve_reaction_kinetics()
    print(f"\n<<< {int(final_answer)} >>>")
