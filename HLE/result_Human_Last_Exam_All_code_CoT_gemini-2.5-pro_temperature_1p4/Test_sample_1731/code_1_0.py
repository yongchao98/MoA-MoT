def find_bose_equilibrium_values():
    """
    This program outlines the derivation and presents the final equilibrium values
    for mean energy and entropy for the Bose case of light quanta (a photon gas),
    as requested by applying principles from large deviation theory.

    The final equations are presented symbolically.
    """
    print("------------------------------------------------------------------")
    print("Equilibrium Energy and Entropy for a Photon Gas (Bose-Einstein Statistics)")
    print("------------------------------------------------------------------")
    print("The equilibrium state is found by maximizing the system's entropy subject to a")
    print("fixed total energy. This procedure, justified by large deviation theory,")
    print("yields the Bose-Einstein distribution. By integrating over all energy states,")
    print("we arrive at the following macroscopic equilibrium values.")
    print("\n--- DERIVED EQUILIBRIUM FORMULAS ---\n")

    # Define symbolic constants for clarity in the output
    E_str = "E"      # Mean Energy
    S_str = "S"      # Entropy
    V_str = "V"      # Volume
    T_str = "T"      # Temperature
    k_B_str = "k_B"  # Boltzmann constant
    hbar_str = "ħ"    # Reduced Planck constant
    c_str = "c"      # Speed of light
    pi_str = "π"     # Pi

    # 1. Equilibrium Mean Energy (E)
    # This is the Stefan-Boltzmann law for the total energy of black-body radiation.
    print("1. Equilibrium Mean Energy (E)")
    print("   The total energy E of a photon gas in a volume V at temperature T is given by:")
    
    numerator_E = f"{pi_str}² * {V_str} * {k_B_str}⁴ * {T_str}⁴"
    denominator_E = f"15 * {hbar_str}³ * {c_str}³"
    
    print(f"\n      {numerator_E}")
    print(f"   {E_str} = {'-' * len(numerator_E)}")
    print(f"      {denominator_E}\n")
    
    print("   The numbers in this equation are:")
    print("   - In the numerator: The exponent '2' on pi (π), and the exponent '4' on the Boltzmann constant (k_B) and temperature (T).")
    print("   - In the denominator: The number '15', and the exponent '3' on the reduced Planck constant (ħ) and the speed of light (c).")
    print("-" * 65)

    # 2. Equilibrium Entropy (S)
    # The entropy is derived from thermodynamic relations (S = (E - F)/T) or directly from the statistical definition.
    print("2. Equilibrium Entropy (S)")
    print("   The total entropy S of the photon gas has a simple relationship with its total energy E and temperature T:")
    
    numerator_S = f"4 * {E_str}"
    denominator_S = f"3 * {T_str}"
    
    print(f"\n      {numerator_S}")
    print(f"   {S_str} = {'-' * len(numerator_S)}")
    print(f"      {denominator_S}\n")
    
    print("   The numbers in this equation are:")
    print("   - In the numerator: The number '4'.")
    print("   - In the denominator: The number '3'.")
    
# Execute the function to display the results
find_bose_equilibrium_values()

# Final answer format for parsing
final_answer_string = "E = (pi^2 * V * k_B^4 * T^4) / (15 * hbar^3 * c^3); S = 4*E / (3*T)"
print(f"\n<<<E = (π² * V * k_B⁴ * T⁴) / (15 * ħ³ * c³), S = 4*E / (3*T)>>>")