import numpy as np
from scipy import constants

def solve_photon_gas_equilibrium():
    """
    Calculates and displays the equilibrium equations for mean energy and
    entropy of a photon gas based on Bose-Einstein statistics and large
    deviation principles.
    """
    print("This script derives the equilibrium values for the mean energy and entropy of a photon gas (light quanta).")
    print("The equilibrium state is the most probable one, found by maximizing entropy for a given mean energy.\n")

    # Define fundamental physical constants
    k_B = constants.k      # Boltzmann constant in J/K
    h = constants.h        # Planck constant in J*s
    c = constants.c        # Speed of light in m/s

    # --- Part 1: Equilibrium Mean Energy <E> ---
    print("--- 1. Equilibrium Mean Energy <E> ---")
    print("The equilibrium mean energy <E> of a photon gas in a volume V at temperature T is given by the Stefan-Boltzmann Law.")
    print("\nThe full theoretical formula is:")
    # Using unicode for pi
    print("<E> = (8 * \u03C0^5 * V * k_B^4 * T^4) / (15 * h^3 * c^3)")

    print("\nAs requested, the numerical constants in this equation are:")
    print(f"  Numerator Coefficient: 8")
    print(f"  Pi exponent: 5")
    print(f"  k_B (Boltzmann constant) exponent: 4")
    print(f"  T (Temperature) exponent: 4")
    print(f"  Denominator Coefficient: 15")
    print(f"  h (Planck constant) exponent: 3")
    print(f"  c (speed of light) exponent: 3")

    # Calculate the Stefan-Boltzmann constant 'a'
    stefan_boltzmann_constant_a = (8 * np.pi**5 * k_B**4) / (15 * h**3 * c**3)

    print("\nThis formula can be simplified to <E>/V = a * T^4, where 'a' is the radiation constant.")
    print("The value of the radiation constant 'a' is:")
    print(f"a = (8 * \u03C0^5 * k_B^4) / (15 * h^3 * c^3) = {stefan_boltzmann_constant_a:.5e} J m^-3 K^-4")

    # --- Part 2: Equilibrium Entropy S ---
    print("\n\n--- 2. Equilibrium Entropy S ---")
    print("The equilibrium entropy S is derived from the thermodynamic relation S = (4/3) * <E> / T.")
    print("\nSubstituting the expression for <E>, the full theoretical formula for entropy S is:")
    print("S = (32 * \u03C0^5 * V * k_B^4 * T^3) / (45 * h^3 * c^3)")

    print("\nAs requested, the numerical constants in this equation are:")
    print(f"  Numerator Coefficient: 32 (= 4/3 * 8)")
    print(f"  Pi exponent: 5")
    print(f"  k_B (Boltzmann constant) exponent: 4")
    print(f"  T (Temperature) exponent: 3")
    print(f"  Denominator Coefficient: 45 (= 3 * 15)")
    print(f"  h (Planck constant) exponent: 3")
    print(f"  c (speed of light) exponent: 3")

    # Calculate the entropy coefficient
    entropy_coefficient_b = (4.0 / 3.0) * stefan_boltzmann_constant_a

    print("\nThis formula can be simplified to S/V = (4/3) * a * T^3.")
    print("The value of the entropy coefficient (4/3)*a is:")
    print(f"(4/3)*a = (32 * \u03C0^5 * k_B^4) / (45 * h^3 * c^3) = {entropy_coefficient_b:.5e} J m^-3 K^-4")


if __name__ == '__main__':
    solve_photon_gas_equilibrium()
    # The final answer is the derivation and explanation provided by the script's output.
    # No single numerical value can be returned without specifying V and T.
    # Thus, the "answer" is the set of derived equations.
    # For a placeholder, let's indicate the derived relationship between S and E.
    final_answer = "S = (4/3) * <E> / T"
    print(f"\n<<<The key relationship found is: {final_answer}>>>")