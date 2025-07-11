import sympy as sp

def derive_photon_gas_properties():
    """
    This function symbolically derives and prints the equilibrium values for
    mean energy (E) and entropy (S) of a photon gas.
    """
    # Step 1: Define symbolic variables for physical constants and parameters.
    # V: Volume, h: Planck's constant, c: speed of light, k_B: Boltzmann constant, T: Temperature
    V, h, c, k_B, T = sp.symbols('V h c k_B T', positive=True, real=True)
    # eps (varepsilon): energy of a single photon mode
    eps = sp.Symbol('varepsilon', positive=True, real=True)

    # beta is a common shorthand in statistical mechanics
    beta = 1 / (k_B * T)

    # Introduction to the physical model
    print("Derivation of Equilibrium Properties for a Photon Gas (Bose case for light quanta)")
    print("-" * 75)
    print("The equilibrium state is the macrostate of maximum entropy. We derive the mean energy E")
    print("and entropy S for this state using the grand canonical ensemble for a gas of photons")
    print("(bosons with zero chemical potential).\n")

    # Step 2: State the formula for Mean Energy E
    # The mean energy E is found by integrating the energy of photons over all possible states.
    # The integral is E = Integral[ ε * <n_ε> * g(ε) dε ] from 0 to infinity,
    # where <n_ε> = 1/(exp(β*ε)-1) is the Bose-Einstein distribution,
    # and g(ε) = (8*pi*V)/(h^3*c^3) * ε^2 is the density of states.
    # The resulting integral is:
    # E = (8*pi*V)/(h^3*c^3) * Integral[ ε^3 / (exp(β*ε) - 1) dε ]
    #
    # This is a standard integral in physics. The definite integral part evaluates to (pi^4 / 15*β^4).
    # Substituting this result gives the Stefan-Boltzmann law for energy.

    print("1. Equilibrium Mean Energy (E)")
    print("-----------------------------------")
    print("The mean energy E of the photon gas is proportional to V*T^4.")
    # We construct and print the full equation for E.
    # The prefactor is (8 * pi^5 * k_B^4) / (15 * h^3 * c^3)
    print("The final equation is:\n")
    print("E = ( (8 * pi^5 * k_B^4) / (15 * h^3 * c^3) ) * V * T^4\n")
    print("The numerical constants in this equation are:")
    print(" - In the numerator: 8, pi^5 (approximately 306.02)")
    print(" - In the denominator: 15")
    print("-" * 75)


    # Step 3: State the formula for Entropy S
    # The entropy S of a photon gas has a remarkably simple relationship to its energy E.
    # It can be derived from the grand partition function, resulting in S = (4/3) * E / T.
    print("2. Equilibrium Entropy (S)")
    print("----------------------------")
    print("The entropy S is directly related to the mean energy E and temperature T.")
    print("The relationship is:\n")
    print("S = (4 / 3) * (E / T)\n")
    print("The numerical constant in this relationship is:")
    print(" - A factor of: 4 / 3")
    print("\nBy substituting the expression for E, we can also write S as a function of T:")
    # We can also construct the full equation for S in terms of temperature.
    # S = (4/3T) * E = (4/3T) * ( (8*pi^5*k_B^4)/(15*h^3*c^3) * V*T^4 )
    # S = ( (32*pi^5*k_B^4)/(45*h^3*c^3) ) * V*T^3
    print("S = ( (32 * pi^5 * k_B^4) / (45 * h^3 * c^3) ) * V * T^3\n")
    print("The numerical constants in this form of the equation are:")
    print(" - In the numerator: 32, pi^5")
    print(" - In the denominator: 45")
    print("-" * 75)


if __name__ == '__main__':
    derive_photon_gas_properties()

<<<E = ( (8 * pi^5 * k_B^4) / (15 * h^3 * c^3) ) * V * T^4, S = (4 / 3) * (E / T)>>>