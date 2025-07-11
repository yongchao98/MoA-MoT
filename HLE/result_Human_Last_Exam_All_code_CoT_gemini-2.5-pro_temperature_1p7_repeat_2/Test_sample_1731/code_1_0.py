import math

def display_equilibrium_values():
    """
    Derives and displays the equilibrium values for mean energy (U) and
    entropy (S) for a photon gas based on Bose-Einstein statistics and
    large deviation principles.
    """
    
    print("This script provides the equilibrium values for mean energy (U) and entropy (S) of a photon gas.")
    print("The derivation follows from maximizing the entropy for a system of bosons with zero chemical potential,")
    print("a principle justified by the Boltzmann-Sanov theorem of large deviation theory.\n")

    print("Step 1: The Bose-Einstein distribution for photons")
    print("Maximizing entropy under the constraint of constant energy leads to the Planck distribution")
    print("for the mean occupation number <n(ε)> of a state with energy ε at temperature T:")
    print("  <n(ε)> = 1 / (exp(ε / (k*T)) - 1)\n")

    print("Step 2: Equilibrium Mean Energy (U)")
    print("The total energy U is found by integrating ε*<n(ε)> over the density of states g(ε) = (8*π*V / (h^3*c^3)) * ε^2.")
    print("  U = Integral from 0 to infinity of [ ε * <n(ε)> * g(ε) ] dε")
    print("Solving this integral leads to the Stefan-Boltzmann law for the total energy:")
    
    # Equation for U
    energy_numerator = "8 * π^5 * V * k^4 * T^4"
    energy_denominator = "15 * h^3 * c^3"
    print("\n--- Equilibrium Mean Energy (U) ---")
    print(f"U = ({energy_numerator}) / ({energy_denominator})\n")

    print("Step 3: Equilibrium Entropy (S)")
    print("Using the thermodynamic relation S = (U - F)/T, where F is the Helmholtz free energy,")
    print("it can be shown for a photon gas that S = 4*U / (3*T).")
    print("Substituting the expression for U gives the equilibrium entropy:")

    # Equation for S
    entropy_numerator = "32 * π^5 * V * k^4 * T^3"
    entropy_denominator = "45 * h^3 * c^3"
    print("\n--- Equilibrium Entropy (S) ---")
    print(f"S = ({entropy_numerator}) / ({entropy_denominator})\n")

    print("Where the constants are:")
    print("  V: Volume of the container")
    print("  T: Absolute temperature")
    print("  k: Boltzmann constant")
    print("  h: Planck constant")
    print("  c: Speed of light in a vacuum")
    print("  π: Mathematical constant Pi")
    
if __name__ == '__main__':
    display_equilibrium_values()
