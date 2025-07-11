import math

def calculate_equilibrium_values():
    """
    Calculates and prints the equilibrium equations for mean energy and entropy
    of a photon gas based on statistical mechanics.
    """
    # Define physical constants in SI units
    h = 6.62607015e-34  # Planck's constant (J*s)
    c = 299792458       # Speed of light (m/s)
    k_B = 1.380649e-23     # Boltzmann constant (J/K)
    PI = math.pi

    # --- Mean Energy Calculation ---
    # The mean energy density u = E/V is given by:
    # u = (8 * pi^5 * k_B^4 / (15 * h^3 * c^3)) * T^4
    energy_coefficient = (8 * PI**5 * k_B**4) / (15 * h**3 * c**3)

    # --- Entropy Calculation ---
    # The entropy density s = S/V is given by:
    # s = (32 * pi^5 * k_B^4 / (45 * h^3 * c^3)) * T^3
    # This can also be expressed as s = (4/3) * u / T
    entropy_coefficient = (32 * PI**5 * k_B**4) / (45 * h**3 * c**3)

    # --- Print the Results ---
    print("The equilibrium values for a gas of light quanta (photons) are functions of temperature T.")
    
    print("\n" + "="*40)
    print("EQUILIBRIUM MEAN ENERGY")
    print("="*40)
    print("The mean energy density, u = E/V (Energy per unit volume), is given by the Stefan-Boltzmann law for black-body radiation.")
    print("\nThe general formula is:")
    print("  u(T) = C_E * T^4")
    print("\nwhere the constant C_E is built from fundamental constants:")
    print("  C_E = (8 * π^5 * k_B^4) / (15 * h^3 * c^3)")
    print(f"\nCalculating the value of this constant:")
    print(f"  π (pi)       = {PI}")
    print(f"  k_B (J/K)    = {k_B}")
    print(f"  h (J*s)      = {h}")
    print(f"  c (m/s)      = {c}")
    print(f"\n  C_E ≈ {energy_coefficient:.6e} J/(m^3 * K^4)")

    print("\nTherefore, the final equation for equilibrium mean energy density is:")
    print(f"  u(T) = {energy_coefficient:.6e} * T^4")

    print("\n" + "="*40)
    print("EQUILIBRIUM ENTROPY")
    print("="*40)
    print("The entropy density, s = S/V (Entropy per unit volume), is derived from the energy.")
    print("\nThe general formula is:")
    print("  s(T) = C_S * T^3")
    print("\nwhere the constant C_S is also built from fundamental constants:")
    print("  C_S = (32 * π^5 * k_B^4) / (45 * h^3 * c^3)")
    print(f"\nCalculating the value of this constant:")
    print(f"  C_S ≈ {entropy_coefficient:.6e} J/(m^3 * K^4)")

    print("\nTherefore, the final equation for equilibrium entropy density is:")
    print(f"  s(T) = {entropy_coefficient:.6e} * T^3")
    
    print("\n" + "="*40)
    print("RELATIONSHIP")
    print("="*40)
    print("Note that energy density and entropy density are related by: s = (4/3) * u / T")

if __name__ == '__main__':
    calculate_equilibrium_values()