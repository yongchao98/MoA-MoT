import math

def calculate_potential_energy_curves():
    """
    Calculates and prints data points for the correct and incorrect
    potential energy curves of the H2+ molecule.

    This script does not run a real quantum chemistry calculation but models
    the potential energy curves to illustrate the concept.
    """

    # Parameters for H2+ in atomic units (Hartree for energy, Bohr for distance)
    De = 0.1026  # Dissociation energy
    Re = 2.00   # Equilibrium bond distance
    a = 0.70    # Morse potential parameter
    E_H = -0.5  # Energy of a separated Hydrogen atom

    # The energy of the dissociated products (H + H+) is E_H
    E_dissociation = E_H

    print("This script models the potential energy surface (PES) for H2+.")
    print("It compares a physically correct curve with an incorrect one that shows")
    print("the unphysical energy drop at long distances you described.\n")
    print(f"{'R (Bohr)':<12} | {'E_correct (Hartree)':<22} | {'E_incorrect (Hartree)':<22}")
    print("-" * 62)

    # Loop through different internuclear distances (R)
    for i in range(10, 101):
        R = i / 10.0

        # 1. Calculate the CORRECT potential energy using the Morse potential
        # E(R) = E_dissociation + De * ((1 - exp(-a(R-Re)))^2 - 1)
        morse_term = (1 - math.exp(-a * (R - Re)))**2 - 1
        E_correct = E_dissociation + De * morse_term

        # 2. Calculate the INCORRECT potential energy
        # This simulates the artifact by adding an unphysical dip at large R.
        # The dip is modeled here with a Gaussian function centered at 6.0 Bohr.
        # The depth of the dip is chosen to go below the equilibrium energy.
        E_incorrect = E_correct
        if R > 3.0:
            dip_magnitude = 0.115
            dip_center = 6.0
            dip_width = 4.0
            dip_term = dip_magnitude * math.exp(-((R - dip_center)**2) / dip_width)
            E_incorrect -= dip_term

        print(f"{R:<12.1f} | {E_correct:<22.6f} | {E_incorrect:<22.6f}")

    # Print the energy at equilibrium for reference
    E_equilibrium = E_dissociation - De
    print("-" * 62)
    print(f"Note: The correct energy at equilibrium (R={Re:.2f} Bohr) is {E_equilibrium:.6f} Hartree.")
    print("The incorrect curve unphysically drops below this value at large R.")


if __name__ == '__main__':
    calculate_potential_energy_curves()
