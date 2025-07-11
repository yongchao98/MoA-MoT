import math

def analyze_demon_apparatus():
    """
    Analyzes the physical conditions for a Maxwell's Demon-like apparatus.

    This function explains why a temperature difference is the crucial parameter
    for a one-way door (a Brownian ratchet) to trap gas on one side. It does
    so by comparing the average kinetic energy of gas molecules in two scenarios:
    one with uniform temperature and one with a temperature difference.
    """

    # Boltzmann constant in Joules per Kelvin
    k_boltzmann = 1.380649e-23

    print("Analysis of Maxwell's Demon with a Brownian Ratchet (One-Way Door)")
    print("="*65)
    print("The one-way door can only cause a net flow of gas if there is an energy")
    print("gradient. This energy is provided by the kinetic energy of gas molecules,")
    print("which is directly proportional to temperature.")
    print("\nLet's calculate the average molecular kinetic energy (KE) using the equation:")
    print("KE = (3/2) * k * T, where k is the Boltzmann constant and T is temperature.")
    print("-" * 65)

    # --- Scenario 1: No Temperature Difference (Thermal Equilibrium) ---
    print("\nSCENARIO 1: UNIFORM TEMPERATURE")
    T1 = 300.0  # Temperature of Chamber 1 in Kelvin
    T2 = 300.0  # Temperature of Chamber 2 in Kelvin
    
    ke1 = (3/2) * k_boltzmann * T1
    ke2 = (3/2) * k_boltzmann * T2

    print(f"If T1 = T2 = {T1} K, the energies are:")
    print(f"  Chamber 1 KE = (3/2) * {k_boltzmann:.4e} J/K * {T1} K = {ke1:.4e} J")
    print(f"  Chamber 2 KE = (3/2) * {k_boltzmann:.4e} J/K * {T2} K = {ke2:.4e} J")
    print("Result: At uniform temperature, molecules on both sides have the same average energy.")
    print("        The door itself jiggles with the same thermal energy. No net flow occurs.")
    print("-" * 65)

    # --- Scenario 2: Temperature Difference ---
    print("\nSCENARIO 2: TEMPERATURE DIFFERENCE")
    T_hot = 400.0  # Temperature of the hot chamber in Kelvin
    T_cold = 300.0 # Temperature of the cold chamber in Kelvin

    ke_hot = (3/2) * k_boltzmann * T_hot
    ke_cold = (3/2) * k_boltzmann * T_cold

    print(f"If T1 = {T_hot} K and T2 = {T_cold} K, the energies are:")
    print(f"  Hot Chamber KE = (3/2) * {k_boltzmann:.4e} J/K * {T_hot} K = {ke_hot:.4e} J")
    print(f"  Cold Chamber KE = (3/2) * {k_boltzmann:.4e} J/K * {T_cold} K = {ke_cold:.4e} J")
    print("Result: Molecules in the hot chamber are more energetic. This energy difference")
    print("        allows the door to be pushed open from the hot side more easily than")
    print("        from the cold side, creating a net flow and trapping gas.")

    print("\nConclusion: The essential experimental parameter is TEMPERATURE.")

if __name__ == '__main__':
    analyze_demon_apparatus()
