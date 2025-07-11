import math

def display_photon_gas_equilibrium_values():
    """
    This function outlines and displays the equilibrium values for mean energy and entropy
    of a photon gas, based on Bose-Einstein statistics and principles of large deviations.
    """
    print("--- Equilibrium Values for a Photon Gas (Bose case) ---")
    print("The equilibrium state is the most probable macrostate, found by maximizing entropy.")
    print("This leads to the Bose-Einstein distribution for photons, which is then used")
    print("to calculate macroscopic properties like mean energy (E) and entropy (S).\n")

    # --- Equilibrium Mean Energy ---
    # The integral for E evaluates to E = (8 * π^5 * k_B^4 / (15 * h^3 * c^3)) * V * T^4
    print("1. Equilibrium Mean Energy (E)")
    print("---------------------------------")
    print("The formula is E = C_E * V * T^4")
    print("where V is volume, T is temperature, and C_E is a constant derived from:")
    print("  C_E = (Numerator) / (Denominator)")
    print("\nBreaking down the constants in the formula for E:")
    
    num_8 = 8
    den_15 = 15
    
    print(f"  Integer in Numerator = {num_8}")
    print( "  Term in Numerator = pi^5 (pi to the power of 5)")
    print( "  Physical constant in Numerator = k_B^4 (Boltzmann's constant to the power of 4)")
    
    print(f"  Integer in Denominator = {den_15}")
    print( "  Physical constant in Denominator = h^3 (Planck's constant to the power of 3)")
    print( "  Physical constant in Denominator = c^3 (speed of light to the power of 3)\n")

    # --- Equilibrium Entropy ---
    # The integral for S evaluates to S = (32 * π^5 * k_B^4 / (45 * h^3 * c^3)) * V * T^3
    print("2. Equilibrium Entropy (S)")
    print("--------------------------")
    print("The formula is S = C_S * V * T^3")
    print("where C_S is a constant derived from:")
    print("  C_S = (Numerator) / (Denominator)")
    print("\nBreaking down the constants in the formula for S:")
    
    num_32 = 32
    den_45 = 45
    
    print(f"  Integer in Numerator = {num_32}")
    print( "  Term in Numerator = pi^5 (pi to the power of 5)")
    print( "  Physical constant in Numerator = k_B^4 (Boltzmann's constant to the power of 4)")
    
    print(f"  Integer in Denominator = {den_45}")
    print( "  Physical constant in Denominator = h^3 (Planck's constant to the power of 3)")
    print( "  Physical constant in Denominator = c^3 (speed of light to the power of 3)\n")

    # --- Relationship between E and S ---
    print("3. Relationship between E and S")
    print("---------------------------------")
    print("A direct and fundamental relationship exists between the equilibrium energy and entropy:")
    
    num_4 = 4
    num_3 = 3
    
    print(f"S = ({num_4} / {num_3}) * E / T")


if __name__ == "__main__":
    display_photon_gas_equilibrium_values()
