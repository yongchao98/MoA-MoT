def solve_bose_light_quanta():
    """
    Calculates and prints the equilibrium values (formulas) for mean energy and entropy
    of a photon gas (Bose case for light quanta).

    The derivation uses the grand canonical partition function for a gas of bosons with
    zero chemical potential, which is consistent with the principles of large deviation theory.

    - Mean Energy <E> is derived from -d(log(Ξ))/dβ.
    - Entropy S is derived from S = (<E> - Ω)/T, where Ω = -k_B*T*log(Ξ).
    """

    # Define symbolic representations of the physical constants and variables.
    k_B = "k_B"  # Boltzmann constant
    h = "h"      # Planck's constant
    c = "c"      # Speed of light
    V = "V"      # Volume
    T = "T"      # Temperature
    pi = "pi"    # Mathematical constant pi

    # The integral ∫[0, inf] x^3 / (e^x - 1) dx evaluates to pi^4 / 15.
    # The density of states g(ε) for photons is proportional to ε^2.
    # These facts are used in the derivation.

    # Assemble the expression for the radiation constant 'a'
    # a = (8 * pi^5 * k_B^4) / (15 * h^3 * c^3)
    energy_coefficient_numerator = f"8 * {pi}^5 * {k_B}^4"
    energy_coefficient_denominator = f"15 * {h}^3 * {c}^3"
    
    # Final formula for Mean Energy <E> (Stefan-Boltzmann Law for energy)
    # The numbers in the equation are 8, 5, 4, 15, 3, 3.
    mean_energy_eq = f"<E> = ({energy_coefficient_numerator} / ({energy_coefficient_denominator})) * {V} * {T}^4"
    
    # The entropy S is related to energy by S = (4/3) * <E> / T
    # So, the coefficient for entropy is (4/3) times the coefficient for <E>/T.
    # (4/3) * (8 * pi^5 * k_B^4 / (15 * h^3 * c^3)) = (32 * pi^5 * k_B^4 / (45 * h^3 * c^3))
    entropy_coefficient_numerator = f"32 * {pi}^5 * {k_B}^4"
    entropy_coefficient_denominator = f"45 * {h}^3 * {c}^3"
    
    # Final formula for Entropy S
    # The numbers in the equation are 32, 5, 4, 45, 3, 3.
    entropy_eq = f"S = ({entropy_coefficient_numerator} / ({entropy_coefficient_denominator})) * {V} * {T}^3"

    print("The equilibrium values are functions of temperature (T) and volume (V).")
    print("They are derived from the partition function, a method whose validity is supported by large deviation theorems.\n")
    print("--- Equilibrium Mean Energy <E> ---")
    print(mean_energy_eq)
    print("\n--- Equilibrium Entropy S ---")
    print(entropy_eq)


if __name__ == "__main__":
    solve_bose_light_quanta()