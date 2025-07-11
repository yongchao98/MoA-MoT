import math

def calculate_rate_constant():
    """
    Calculates the reaction rate constant based on the Energetic Span Model
    from the provided reaction energy profile.
    """
    # Step 1: Define constants and input values from the diagram
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34   # Planck constant in J*s
    R = 8.314462         # Ideal gas constant in J/(mol*K)
    T = 298              # Temperature in Kelvin

    # Energies from the provided diagram in kJ/mol
    G_TS_eli = 18.0
    G_Cat_GB = -61.0

    # Step 2: Identify the rate-determining step and calculate the energetic span (delta_E)
    print("--- Identifying the Rate-Determining Step and Energetic Span ---")
    print("According to the Energetic Span Model, the rate is determined by the highest energy transition state (TDTS) and the lowest energy intermediate preceding it.")
    print(f"The highest energy transition state (TDTS) is TS-eli at G = {G_TS_eli} kJ/mol.")
    print(f"The lowest energy intermediate preceding TS-eli is Cat-GB at G = {G_Cat_GB} kJ/mol.")
    print("Therefore, the rate-determining step is the elimination step associated with TS-eli (P5 -> P).")
    
    delta_E_kJ = G_TS_eli - G_Cat_GB
    print("\nThe energetic span (δE) is calculated as:")
    print(f"δE = G(TS-eli) - G(Cat-GB) = {G_TS_eli} kJ/mol - ({G_Cat_GB} kJ/mol) = {delta_E_kJ} kJ/mol")

    # Convert delta_E to J/mol for calculation
    delta_E_J = delta_E_kJ * 1000

    # Step 3: Calculate the rate constant (k) using the Eyring equation
    # k = (k_B * T / h) * exp(-delta_E / (R * T))
    pre_exponential_factor = (k_B * T) / h
    RT = R * T
    exponential_term = math.exp(-delta_E_J / RT)
    
    k_per_second = pre_exponential_factor * exponential_term

    # Step 4: Convert k to hours^-1 and round to two significant figures
    k_per_hour = k_per_second * 3600
    
    # Rounding to two significant figures
    k_final = float(f"{k_per_hour:.1e}")

    print("\n--- Calculating the Reaction Rate Constant (k) ---")
    print("The rate constant is calculated using the Eyring equation: k = (k_B*T/h) * exp(-δE/RT)")
    print(f"Using δE = {delta_E_J} J/mol, T = {T} K:")
    print(f"k = {k_per_second:.4f} s^-1")
    print(f"Converting to hours^-1: k = {k_per_second:.4f} * 3600 = {k_per_hour:.2f} hours^-1")
    print(f"\nThe calculated reaction rate constant to two significant figures is: {k_final} hours^-1")
    
    # Final answer block
    # Convert scientific notation to integer for final format if possible
    final_answer = int(k_final) if k_final.is_integer() else k_final
    
    # Since 2.9e2 is 290, we can represent it as an integer.
    if 'e' in str(k_final):
        base, exponent = str(k_final).split('e')
        final_answer = int(float(base) * (10**int(exponent)))
        
    print(f"\n<<<k = {final_answer}>>>")

calculate_rate_constant()