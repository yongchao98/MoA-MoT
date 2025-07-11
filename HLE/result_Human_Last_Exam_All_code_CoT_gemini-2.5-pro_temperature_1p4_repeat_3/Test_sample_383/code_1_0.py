import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_pes(distances):
    """
    Calculates the potential energy surface of H2+ using Restricted Hartree-Fock (RHF).

    For a one-electron system like H2+, RHF is an exact method (within the basis set)
    and does not suffer from the self-interaction error that causes the energy to
    drop incorrectly at long distances in many DFT functionals.

    The expected dissociation energy is -0.5 Hartree, the energy of an isolated
    hydrogen atom.
    """
    print("Calculating H2+ Potential Energy Curve using Hartree-Fock")
    print("---------------------------------------------------------")
    print("Distance (Angstrom)    Total Energy (Hartree)")
    print("---------------------------------------------------------")

    for r in distances:
        # Define the molecule: H2+
        # Set charge=1 and spin=1 (one unpaired electron -> doublet state)
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, r)]
        ]
        mol.charge = 1
        mol.spin = 1
        mol.basis = 'cc-pvdz'
        mol.build()

        # Perform a Restricted Hartree-Fock (RHF) calculation
        mf = scf.RHF(mol).run(quiet=True)

        # Print the distance and calculated total energy
        # The dissociation limit is H + H+, which has an energy of a single H atom (-0.5 Hartree)
        # Note: total_energy = electronic_energy + nuclear_repulsion_energy
        total_energy = mf.e_tot
        print(f"{r:<20.4f} {total_energy:<20.10f}")

    print("---------------------------------------------------------")
    print("Note: As distance increases, the energy correctly approaches -0.5 Hartree.")


if __name__ == '__main__':
    # Define a range of bond distances to scan, from equilibrium to dissociation
    # Equilibrium distance is around 1.06 Angstroms
    bond_distances = np.concatenate([
        np.linspace(0.5, 2.0, 16),
        np.linspace(2.2, 5.0, 15),
        np.array([6.0, 7.0, 8.0, 10.0, 20.0])
    ])
    
    calculate_h2_plus_pes(bond_distances)
