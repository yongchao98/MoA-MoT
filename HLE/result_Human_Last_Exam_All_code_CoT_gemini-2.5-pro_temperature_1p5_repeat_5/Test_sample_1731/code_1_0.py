import math

def calculate_photon_gas_properties():
    """
    Calculates and displays the equilibrium formulas for mean energy and entropy
    of a photon gas based on statistical mechanics.
    """

    # We use symbolic placeholders for physical constants and variables,
    # as no numerical values for Volume (V) or Temperature (T) are given.
    k_B = "k_B"  # Boltzmann constant
    h = "h"      # Planck's constant
    c = "c"      # Speed of light
    V = "V"      # Volume
    T = "T"      # Temperature

    print("--- Equilibrium Properties of a Photon Gas (Bose Case) ---")
    
    # --- Mean Energy <E> ---
    print("\n1. Equilibrium Mean Energy <E>")
    print("-" * 35)

    # The calculation involves a standard definite integral:
    # integral from 0 to infinity of [x^3 / (exp(x) - 1)] dx = pi^4 / 15
    integral_val_energy = math.pi**4 / 15

    # The full formula for <E> is:
    # <E> = (8 * pi * V * (k_B*T)^4 / (h^3 * c^3)) * (pi^4 / 15)
    # <E> = (8 * pi^5 * V * k_B^4 / (15 * h^3 * c^3)) * T^4
    
    # Calculate the numerical part of the coefficient
    energy_coeff_numeric = (8 * math.pi**5) / 15
    
    print("The equilibrium mean energy <E> is given by the Stefan-Boltzmann law for total energy.")
    print("\nThe equation for <E> is:")
    print(f"<E> = C_E * T^4")
    
    print("\nWhere the energy coefficient C_E is:")
    # Print each part of the coefficient clearly
    print(f"C_E = (8 * π^5 * V * {k_B}^4) / (15 * {h}^3 * {c}^3)")
    print(f"The numerical part of the constant is (8 * π^5 / 15) ≈ {energy_coeff_numeric:.4f}")

    print("\nFinal Equation for Mean Energy:")
    print(f"      8 * π^5 * V * {k_B}^4")
    print(f"<E> = ------------------- * T^4")
    print(f"        15 * {h}^3 * {c}^3")


    # --- Entropy S ---
    print("\n\n2. Equilibrium Entropy S")
    print("-" * 35)
    
    # For a photon gas, there's a simple relation: S = (4/3) * <E> / T
    # S = (4/3) * (1/T) * [ (8*pi^5*V*k_B^4)/(15*h^3*c^3) * T^4 ]
    # S = (32 * pi^5 * V * k_B^4 / (45 * h^3 * c^3)) * T^3
    
    # Calculate the numerical part of the coefficient
    entropy_coeff_numeric = (32 * math.pi**5) / 45
    
    print("The equilibrium entropy S is derived from the mean energy via the relation S = 4*<E>/(3*T).")
    print("\nThe equation for S is:")
    print(f"S = C_S * T^3")

    print("\nWhere the entropy coefficient C_S is:")
    print(f"C_S = (32 * π^5 * V * {k_B}^4) / (45 * {h}^3 * {c}^3)")
    print(f"The numerical part of the constant is (32 * π^5 / 45) ≈ {entropy_coeff_numeric:.4f}")

    print("\nFinal Equation for Entropy:")
    print(f"      32 * π^5 * V * {k_B}^4")
    print(f"S =  -------------------- * T^3")
    print(f"        45 * {h}^3 * {c}^3")

if __name__ == '__main__':
    calculate_photon_gas_properties()
