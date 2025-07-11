import math

def calculate_photon_gas_properties():
    """
    Calculates and prints the equilibrium mean energy and entropy for a photon gas
    based on the principles of Bose-Einstein statistics.
    """
    # 1. Define physical constants
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    h = 6.62607015e-34   # Planck constant (J*s)
    c = 299792458        # Speed of light (m/s)
    pi = math.pi

    # 2. Set example parameters for the calculation
    T = 1500.0  # Temperature in Kelvin (e.g., surface of a red dwarf star)
    V = 1.0     # Volume in cubic meters

    # 3. Calculate the constant term for the energy formula.
    # This term is the Stefan-Boltzmann law's energy density constant (sigma').
    # sigma_prime = (8 * pi^5 * k_B^4) / (15 * h^3 * c^3)
    sigma_prime_numerator = 8 * (pi**5) * (k_B**4)
    sigma_prime_denominator = 15 * (h**3) * (c**3)
    sigma_prime = sigma_prime_numerator / sigma_prime_denominator

    # 4. Calculate Equilibrium Mean Energy (E)
    # The formula is E = V * sigma_prime * T^4
    E = V * sigma_prime * (T**4)

    # 5. Calculate Equilibrium Entropy (S)
    # The formula is S = (4/3) * E / T
    S = (4 / 3) * E / T

    # 6. Print the results in a detailed format
    print("Equilibrium Values for a Photon Gas (Bose Case)")
    print("-" * 50)
    print(f"Assuming Temperature T = {T} K and Volume V = {V} m^3\n")

    # --- Mean Energy ---
    print("1. Equilibrium Mean Energy (E)")
    print("Formula: E = (8 * π^5 * k_B^4 / (15 * h^3 * c^3)) * V * T^4")
    print("Breaking down the final equation:")
    print(f"  Constant Term (8 * π^5 * k_B^4 / (15 * h^3 * c^3)) = {sigma_prime:.4e} J/(m^3*K^4)")
    print(f"  Volume (V) = {V} m^3")
    print(f"  Temperature (T) = {T} K")
    print("\nFinal Equation with numbers:")
    print(f"E = {sigma_prime:.4e} * {V} * {T}^4")
    print(f"E = {E:.4e} Joules")
    print("-" * 50)

    # --- Entropy ---
    print("2. Equilibrium Entropy (S)")
    print("Formula: S = (4/3) * E / T")
    print("Breaking down the final equation:")
    print(f"  Factor = {4/3:.4f}")
    print(f"  Mean Energy (E) = {E:.4e} J")
    print(f"  Temperature (T) = {T} K")
    print("\nFinal Equation with numbers:")
    print(f"S = (4/3) * {E:.4e} / {T}")
    print(f"S = {S:.4e} J/K")
    print("-" * 50)

if __name__ == '__main__':
    calculate_photon_gas_properties()