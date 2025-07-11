import numpy as np

def calculate_photon_gas_properties(T, V):
    """
    Calculates the equilibrium mean energy and entropy for a photon gas (Bose case)
    based on large deviation principles and statistical mechanics.

    Args:
        T (float): Temperature in Kelvin (K).
        V (float): Volume in cubic meters (m^3).
    """
    # Define physical constants
    h = 6.62607015e-34  # Planck's constant (J·s)
    c = 2.99792458e8    # Speed of light (m/s)
    k_B = 1.380649e-23   # Boltzmann constant (J/K)

    print(f"Calculating for Temperature T = {T} K and Volume V = {V} m^3\n")

    # --- Mean Energy (U) Calculation ---
    U_numerator = 8 * (np.pi**5) * V * ((k_B * T)**4)
    U_denominator = 15 * (h**3) * (c**3)
    U = U_numerator / U_denominator

    print("--- Equilibrium Mean Energy (U) ---")
    print("Equation: U = (8 * pi^5 * V * (k_B * T)^4) / (15 * h^3 * c^3)")
    print("Plugging in the values:")
    print(f"U = (8 * {np.pi**5:.4f} * {V:.2e} * ({k_B:.4e} * {T:.1e})^4) / (15 * ({h:.4e})^3 * ({c:.4e})^3)")
    print(f"Result: U = {U:.6e} Joules\n")

    # --- Entropy (S) Calculation ---
    S_numerator = 32 * (np.pi**5) * V * k_B * ((k_B * T)**3)
    S_denominator = 45 * (h**3) * (c**3)
    S = S_numerator / S_denominator

    print("--- Equilibrium Entropy (S) ---")
    print("Equation: S = (32 * pi^5 * V * k_B^4 * T^3) / (45 * h^3 * c^3)")
    print("Plugging in the values:")
    print(f"S = (32 * {np.pi**5:.4f} * {V:.2e} * {k_B:.4e} * ({k_B:.4e} * {T:.1e})^3) / (45 * ({h:.4e})^3 * ({c:.4e})^3)")
    print(f"Result: S = {S:.6e} J/K")
    
    # As a check, verify S = 4U / (3T)
    # S_check = (4 * U) / (3 * T)
    # print(f"\nVerification check: 4*U/(3*T) = {S_check:.6e} J/K (matches S)")


# Example Usage: Calculate for a 1 cm^3 cavity at 1500 K (like a furnace)
temperature = 1500.0  # Kelvin
volume = 1e-6       # 1 cm^3 = 1e-6 m^3

calculate_photon_gas_properties(temperature, volume)

final_answer = "The equilibrium values are given by the functions U(T,V) and S(T,V) derived and implemented in the code."
<<<The equilibrium values are given by the functions U(T,V) and S(T,V) derived and implemented in the code.>>>