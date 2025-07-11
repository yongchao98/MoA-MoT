import math
import scipy.constants

def calculate_photon_gas_properties():
    """
    This function calculates the coefficients for the equilibrium energy and
    entropy of a photon gas based on fundamental constants. It then prints
    the final equations in a clear format.
    """

    # --- Fundamental constants from scipy ---
    k = scipy.constants.k      # Boltzmann constant (J/K)
    h = scipy.constants.h      # Planck constant (J*s)
    c = scipy.constants.c      # Speed of light (m/s)
    pi = scipy.constants.pi    # Pi

    # --- Derivation of Coefficients ---
    # The integral ∫ x^3 / (e^x - 1) dx from 0 to ∞ is π^4 / 15.
    # Mean Energy E = ( (8*π*V) / (h^3*c^3) ) * ∫ ε^3 / (exp(ε/(kT)) - 1) dε
    # After substitution x = ε/(kT), this leads to:
    # E = ( (8*π*V) / (h^3*c^3) ) * (kT)^4 * (π^4/15)
    # E = [ (8*π^5*k^4) / (15*h^3*c^3) ] * V * T^4
    energy_coeff = (8 * pi**5 * k**4) / (15 * h**3 * c**3)

    # Entropy S and Energy E for a photon gas are related by S = (4/3) * E / T.
    # S = (4/3) * (energy_coeff * V * T^4) / T
    # S = [ (4/3) * energy_coeff ] * V * T^3
    # S = [ (32*π^5*k^4) / (45*h^3*c^3) ] * V * T^3
    entropy_coeff = (4 / 3) * energy_coeff

    # --- Output the results ---
    print("Based on the principles of large deviations and statistical mechanics for a Bose gas of photons:")
    print("-" * 80)
    
    # Equilibrium Mean Energy
    print("The equilibrium mean energy (E) of a photon gas in a volume V at temperature T is:")
    print(f"E = [ (8 * π^5 * k^4) / (15 * h^3 * c^3) ] * V * T^4")
    # Output the equation with the calculated numeric coefficient
    # Each number is outputted in the formatted string.
    print("\nNumerically, the equation for energy is:")
    final_energy_eq = f"E = {energy_coeff:.4e} * V * T^4  (in Joules)"
    print(final_energy_eq)
    
    print("-" * 80)

    # Equilibrium Entropy
    print("The equilibrium entropy (S) of a photon gas in a volume V at temperature T is:")
    print(f"S = [ (32 * π^5 * k^4) / (45 * h^3 * c^3) ] * V * T^3")
    # Output the equation with the calculated numeric coefficient
    print("\nNumerically, the equation for entropy is:")
    final_entropy_eq = f"S = {entropy_coeff:.4e} * V * T^3  (in Joules/Kelvin)"
    print(final_entropy_eq)
    
    print("-" * 80)


if __name__ == '__main__':
    calculate_photon_gas_properties()
