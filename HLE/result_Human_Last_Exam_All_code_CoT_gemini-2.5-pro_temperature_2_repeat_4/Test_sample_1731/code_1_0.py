import math

def display_photon_gas_equilibrium_values():
    """
    This function derives and displays the equilibrium values for mean energy and entropy
    for the Bose case of light quanta (a photon gas), based on principles of statistical
    mechanics, which are formally justified by large deviation theorems.
    """

    print("Derivation of Equilibrium Values for a Photon Gas")
    print("=" * 50)
    print("The system is a gas of photons (light quanta) in a volume V at a temperature T.")
    print("Photons are bosons with zero chemical potential.")
    print("The equilibrium state is found by maximizing entropy subject to a fixed mean energy.")
    print("\nThis leads to the Bose-Einstein distribution for photons: n(e) = 1 / (exp(e/(k_B*T)) - 1)")
    print("\nIntegrating over the density of states for photons yields the following results.")
    print("-" * 50)

    # --- Variable Definitions for the equations ---
    # These are string representations for the final output.
    energy = "<E>"
    entropy = "S"
    volume = "V"
    temp = "T"
    boltzmann_const = "k_B"
    planck_const = "h"
    speed_of_light = "c"
    pi_symbol = "pi"

    # --- Mean Energy ---
    print("\n1. Equilibrium Mean Energy ({})".format(energy))
    print("\nThe mean energy of the photon gas is given by:")
    print(
        "{energy} = (8 * {pi}^5 * {V} * {k_B}^4 / (15 * {h}^3 * {c}^3)) * {T}^4"
        .format(
            energy=energy,
            pi=pi_symbol,
            V=volume,
            k_B=boltzmann_const,
            h=planck_const,
            c=speed_of_light,
            T=temp
        )
    )
    print("\nThis is the Stefan-Boltzmann law for the total energy of black-body radiation.")


    # --- Entropy ---
    print("\n" + "-" * 50)
    print("\n2. Equilibrium Entropy ({})".format(entropy))
    print("\nThe entropy of the photon gas, found via the thermodynamic relation S = (4/3) * <E> / T, is:")
    print(
        "{S} = (32 * {pi}^5 * {V} * {k_B}^4 / (45 * {h}^3 * {c}^3)) * {T}^3"
        .format(
            S=entropy,
            pi=pi_symbol,
            V=volume,
            k_B=boltzmann_const,
            h=planck_const,
            c=speed_of_light,
            T=temp
        )
    )
    print("\n" + "=" * 50)

# Execute the function to display the results
if __name__ == "__main__":
    display_photon_gas_equilibrium_values()