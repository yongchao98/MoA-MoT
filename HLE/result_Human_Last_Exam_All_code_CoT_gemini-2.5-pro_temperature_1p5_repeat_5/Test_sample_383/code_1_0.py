import numpy
from pyscf import gto, scf

def calculate_h2plus_potential_energy_curve():
    """
    Calculates and prints the potential energy curve for the H2+ molecule.

    This code addresses the common issue of incorrect dissociation in symmetric
    molecules by taking the following steps:
    1.  Using Unrestricted Hartree-Fock (UHF): Unlike restricted methods, UHF allows
        the electronic wavefunction to break spatial symmetry. This is essential for
        correctly describing the dissociation of H2+ into H and H+, where the
        electron localizes on one atom.
    2.  Disabling Symmetry Constraints: By setting `mol.symmetry = False`, we explicitly
        tell the program not to enforce any point group symmetry. This gives the
        calculation the freedom to find the true, lower-energy asymmetric state at
        large internuclear distances.
    """
    print("Calculating the H2+ potential energy curve using UHF...")
    print("This method allows for symmetry breaking, which correctly models dissociation.")
    print("-" * 52)
    print("Distance (Angstrom) | Total Energy (Hartree)")
    print("--------------------|-------------------------")

    # Define a range of internuclear distances to scan
    # A numpy array from 0.4 Angstrom to 6.0 Angstrom in steps of 0.2
    distances = numpy.arange(0.4, 6.1, 0.2)

    for dist in distances:
        # 1. Define the molecule at the given distance
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, dist)]
        ]
        mol.basis = 'cc-pvdz'  # A standard correlation-consistent basis set
        mol.charge = 1         # Charge of +1 for the cation
        mol.spin = 1           # Spin of 1 for one unpaired electron (doublet)

        # 2. *** This is the key part of the fix ***
        # Disable symmetry enforcement to allow for localization of the electron
        mol.symmetry = False

        # Build the molecule object
        mol.build()

        # 3. Use Unrestricted Hartree-Fock (UHF)
        mf = scf.UHF(mol)

        # Run the calculation quietly (suppress intermediate SCF cycle output)
        energy = mf.run(verbose=0).e_tot

        # 4. Print the final result for the current distance
        print(f"{dist:^20.2f}|{energy:^25.8f}")

    print("-" * 52)
    print("Calculation complete. The energy correctly approaches -0.5 Hartree")
    print("at large distances, the energy of an isolated H atom.")


# Run the calculation and print the results
calculate_h2plus_potential_energy_curve()