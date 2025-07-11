import numpy as np

def calculate_photon_gas_properties():
    """
    Calculates and prints the equilibrium mean energy and entropy of a photon gas.
    
    This script applies the results from statistical mechanics for a gas of bosons
    (photons) where the number of particles is not conserved. The underlying principles,
    such as maximizing entropy for a given mean energy, are physical applications of
    large deviation theorems.
    
    The script first presents the derived symbolic formulas and then calculates
    numerical values for an example case.
    """

    # --- 1. Define physical constants ---
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34   # Planck constant in J*s
    c = 2.99792458e8     # Speed of light in m/s
    pi = np.pi

    # --- 2. Define example system parameters ---
    V = 1.0  # Volume in cubic meters (m^3)
    T = 300.0 # Temperature in Kelvin (K)

    # --- 3. Introduction and Explanation ---
    print("This script calculates the equilibrium mean energy and entropy for a photon gas.")
    print("The results are derived by maximizing the Bose-Einstein entropy under the constraint")
    print("of a fixed average energy, a principle underpinned by large deviation theory.")
    print(f"\nSystem Parameters: Volume V = {V:.1f} m^3, Temperature T = {T:.1f} K")
    print("-" * 70)

    # --- 4. Equilibrium Mean Energy Calculation ---
    print("\n1. Equilibrium Mean Energy (<E>):")
    print("\nThe formula for the mean energy of a photon gas in a volume V at temperature T is:")
    print("   <E> = (8 * π^5 * k_B^4 / (15 * h^3 * c^3)) * V * T^4")

    # The constant term is related to the Stefan-Boltzmann law energy density (u = a*T^4)
    energy_const = (8 * pi**5 * k_B**4) / (15 * h**3 * c**3)
    mean_energy = energy_const * V * T**4
    
    print("\n   Plugging in the numerical values:")
    print(f"   <E> = (8 * ({pi:.5f})^5 * ({k_B:.6e})^4 / (15 * ({h:.6e})^3 * ({c:.6e})^3)) * {V:.1f} * ({T:.1f})^4")
    print(f"\n   The calculated mean energy is:")
    print(f"   <E> = {mean_energy:.6e} Joules")
    print("-" * 70)

    # --- 5. Equilibrium Entropy Calculation ---
    print("\n2. Equilibrium Entropy (S):")
    print("\nThe formula for the entropy is derived from the thermodynamic relation S = (4/3) * <E> / T:")
    print("   S = (32 * π^5 * k_B^4 / (45 * h^3 * c^3)) * V * T^3")

    # We can calculate entropy directly from the mean energy
    entropy = (4.0 / 3.0) * mean_energy / T
    entropy_const = (32 * pi**5 * k_B**4) / (45 * h**3 * c**3)

    print("\n   Plugging in the numerical values:")
    print(f"   S = (32 * ({pi:.5f})^5 * ({k_B:.6e})^4 / (45 * ({h:.6e})^3 * ({c:.6e})^3)) * {V:.1f} * ({T:.1f})^3")
    print(f"\n   The calculated entropy is:")
    print(f"   S = {entropy:.6e} J/K")
    print("-" * 70)

if __name__ == '__main__':
    calculate_photon_gas_properties()