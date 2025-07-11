import math

def display_photon_gas_equations():
    """
    Derives and displays the equilibrium values for mean energy (U) and entropy (S)
    for a photon gas (Bose case) based on statistical mechanics, which is underpinned
    by large deviation principles.
    """

    # --- Variable Definitions for the equations ---
    # These are symbolic representations for the final output.
    symbols = {
        'U': 'U', # Mean Energy
        'S': 'S', # Entropy
        'V': 'V', # Volume
        'T': 'T', # Temperature
        'k_B': 'k_B', # Boltzmann constant
        'hbar': 'ħ',  # Reduced Planck constant (h / 2π)
        'c': 'c',   # Speed of light in vacuum
        'pi': 'π'    # Mathematical constant pi
    }

    # --- Mean Energy (U) Derivation ---
    # The derivation involves integrating ε * g(ε) / (exp(ε/(k_B*T)) - 1) from 0 to infinity.
    # The solution to the dimensionless part of the integral, ∫x^3/(e^x-1)dx, is π^4/15.
    energy_numerator = f"{symbols['pi']}^2 * {symbols['V']} * {symbols['k_B']}^4 * {symbols['T']}^4"
    energy_denominator_num = 15
    energy_denominator = f"{energy_denominator_num} * {symbols['hbar']}^3 * {symbols['c']}^3"
    
    print("--- Equilibrium Mean Energy (U) for a Photon Gas ---")
    print("The equilibrium mean energy is derived by integrating the energy contribution of all states,")
    print("weighted by the Bose-Einstein distribution for photons.")
    print("\nThe final equation is:")
    print(f"  {symbols['U']} = ( {energy_numerator} ) / ( {energy_denominator} )")
    
    print(f"\n> The number '{energy_denominator_num}' in the denominator originates from the solution to the standard definite integral")
    print(f"  ∫[x^3 / (e^x - 1)]dx from 0 to ∞, which evaluates to π^4 / {energy_denominator_num}.")

    # --- Entropy (S) Derivation ---
    # Derived from the thermodynamic relation S = (U - F) / T, where for a photon gas F = -U/3.
    # This simplifies to S = 4*U / (3*T).
    # Substituting the expression for U results in 4 / (3 * 15) = 4 / 45.
    entropy_numerator = f"4 * {symbols['pi']}^2 * {symbols['V']} * {symbols['k_B']}^4 * {symbols['T']}^3"
    entropy_denominator_num = 45 # This comes from 3 * 15 from the S = 4U/(3T) relation
    entropy_denominator = f"{entropy_denominator_num} * {symbols['hbar']}^3 * {symbols['c']}^3"

    print("\n\n--- Equilibrium Entropy (S) for a Photon Gas ---")
    print("The equilibrium entropy is found using the thermodynamic relation between energy and entropy.")
    print("For a photon gas, this relation is S = 4*U / (3*T).")
    print("\nThe final equation is:")
    print(f"  {symbols['S']} = ( {entropy_numerator} ) / ( {entropy_denominator} )")

    print(f"\n> The number '{entropy_denominator_num}' in the denominator is the result of the factor 4/3 from the")
    print(f"  thermodynamic relation multiplied by the factor 1/{energy_denominator_num} from the energy equation (3 * {energy_denominator_num} = {entropy_denominator_num}).")

# Execute the function to display the results
if __name__ == "__main__":
    display_photon_gas_equations()
