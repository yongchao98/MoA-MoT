import math

def solve_photon_gas_equilibrium():
    """
    This script explains the derivation and presents the equilibrium values for
    mean energy and entropy of a photon gas (Bose case for light quanta) based
    on statistical mechanics, whose foundational principles are supported by
    large deviation theory.
    """

    print("--- Theoretical Framework ---")
    print("The equilibrium state of the photon gas is the state of maximum entropy.")
    print("This is justified by Large Deviation Theory, which states that deviations")
    print("from this state are exponentially improbable.\n")
    print("For a gas of photons (bosons with non-conserved number), the chemical potential is 0.")
    print("The average number of photons n(e) in an energy state 'e' is given by the Bose-Einstein distribution:")
    print("  <n(e)> / g(e) = 1 / (exp(e / (k_B * T)) - 1)")
    print("where g(e) is the density of states, k_B is the Boltzmann constant, and T is temperature.\n")
    print("The density of states g(e) for photons in a volume V is:")
    print("  g(e)de = (8 * pi * V / (h^3 * c^3)) * e^2 de\n")

    # --- Mean Energy Calculation ---
    print("--- 1. Equilibrium Mean Energy (E) ---")
    print("The mean energy E is found by integrating the energy 'e' multiplied by the number")
    print("of photons n(e) = <n(e)> over all energies:")
    print("  E = Integral from 0 to inf [ e * g(e) * (1 / (exp(e / (k_B * T)) - 1)) de ]")
    print("  E = Integral [ (8 * pi * V / (h^3 * c^3)) * e^3 / (exp(e / (k_B * T)) - 1) de ]\n")
    print("This integral can be solved using the substitution x = e / (k_B * T).")
    print("The integral part reduces to a standard form: Integral [ x^3 / (e^x - 1) dx ]")
    
    integral_val_E = (math.pi**4) / 15
    print(f"The value of this standard integral is pi^4 / 15 approx {integral_val_E:.4f}\n")

    print("After performing the integration, we get the final equation for the mean energy E (Stefan-Boltzmann Law):")
    
    # Constants for the energy equation
    num_E_numerator = 8
    num_E_denominator = 15
    
    print("E = ({} * \u03C0\u2075 * k_B\u2074 * T\u2074 * V) / ({} * h\u00b3 * c\u00b3)".format(num_E_numerator, num_E_denominator))
    print("Each number in the final equation is:")
    print(f"Numerator coefficient: {num_E_numerator}")
    print(f"Denominator coefficient: {num_E_denominator}\n")
    
    # --- Entropy Calculation ---
    print("--- 2. Equilibrium Entropy (S) ---")
    print("The entropy S can be derived from the thermodynamic relation TS = E + PV.")
    print("For a photon gas, the pressure P relates to energy density E/V by P = E / (3V).")
    print("Therefore, TS = E + (E/3) = 4E/3, which means S = 4E / (3T).\n")
    print("Substituting the derived expression for E gives:")
    
    # Constants for the entropy equation
    # S = (4 / (3T)) * (8 * pi^5 * k_B^4 * T^4 * V) / (15 * h^3 * c^3)
    # S = (32 * pi^5 * k_B^4 * T^3 * V) / (45 * h^3 * c^3)
    num_S_numerator = 32
    num_S_denominator = 45

    print("S = ({} * \u03C0\u2075 * k_B\u2074 * T\u00b3 * V) / ({} * h\u00b3 * c\u00b3)".format(num_S_numerator, num_S_denominator))
    print("Each number in the final equation is:")
    print(f"Numerator coefficient: {num_S_numerator}")
    print(f"Denominator coefficient: {num_S_denominator}\n")

if __name__ == '__main__':
    solve_photon_gas_equilibrium()