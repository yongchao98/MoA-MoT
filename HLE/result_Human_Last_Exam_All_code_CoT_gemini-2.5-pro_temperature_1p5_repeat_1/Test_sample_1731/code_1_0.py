def display_equilibrium_values():
    """
    Derives and displays the equilibrium values for mean energy and entropy
    for a photon gas (Bose case of light quanta) based on the principles
    of statistical mechanics, which are justified by large deviation theory.

    The derivation proceeds as follows:
    1. The equilibrium state is the one that maximizes the Boltzmann entropy
       S = k_B * ln(W) for bosons, under the constraint of a fixed total energy E.
    2. This maximization yields the Bose-Einstein distribution for the mean
       occupation number <n_j> of an energy state ε_j:
       <n_j> = 1 / (exp(ε_j / (k_B*T)) - 1)
    3. The total mean energy (E) is calculated by integrating the energy of
       all states, weighted by this distribution, over the density of states
       g(ε) for photons in a 3D volume V:
       E = Integral[ε * <n(ε)> * g(ε) dε] from 0 to infinity.
       The density of states is g(ε) = (8 * π * V / (h^3 * c^3)) * ε^2.
    4. Performing this integration leads to the formula for E shown below, which is
       a form of the Stefan-Boltzmann law.
       The key integral is Integral[x^3 / (e^x - 1) dx] = π^4 / 15.
    5. The entropy (S) is found using the thermodynamic relation for a photon gas,
       S = (4/3) * (E/T). Substituting the expression for E gives the formula for S.
    """

    # Define the components of the final equations as strings
    mean_energy_numerator = "8 * π^5 * V * k_B^4 * T^4"
    mean_energy_denominator = "15 * h^3 * c^3"

    entropy_numerator = "32 * π^5 * V * k_B^4 * T^3"
    entropy_denominator = "45 * h^3 * c^3"

    print("--- Equilibrium Values for Photon Gas ---")
    print("\nPhysical constants used in the equations:")
    print("  k_B: Boltzmann constant")
    print("  h: Planck's constant")
    print("  c: Speed of light")
    print("  π: Pi")
    print("\nSystem variables:")
    print("  T: Absolute temperature")
    print("  V: Volume of the system")
    
    print("\n--- Equilibrium Mean Energy (E) ---")
    print("The equation for the mean energy is:")
    print(f"E = ({mean_energy_numerator}) / ({mean_energy_denominator})")
    
    print("\nDecomposing the mean energy equation:")
    print(f"Numerator: {mean_energy_numerator}")
    print("It includes the numerical factors: 8, 5 (as π^5), 4 (as T^4 and k_B^4)")
    print(f"Denominator: {mean_energy_denominator}")
    print("It includes the numerical factors: 15, 3 (as h^3 and c^3)")

    print("\n--- Equilibrium Entropy (S) ---")
    print("The equation for the entropy is:")
    print(f"S = ({entropy_numerator}) / ({entropy_denominator})")

    print("\nDecomposing the entropy equation:")
    print(f"Numerator: {entropy_numerator}")
    print("It includes the numerical factors: 32, 5 (as π^5), 4 (as k_B^4), 3 (as T^3)")
    print(f"Denominator: {entropy_denominator}")
    print("It includes the numerical factors: 45, 3 (as h^3 and c^3)")


# Execute the function to display the results
display_equilibrium_values()
<<<
E = (8 * π^5 * V * k_B^4 * T^4) / (15 * h^3 * c^3)
S = (32 * π^5 * V * k_B^4 * T^3) / (45 * h^3 * c^3)
>>>