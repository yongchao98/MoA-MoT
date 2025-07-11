import numpy as np
from pyscf import gto, scf, dft

def calculate_h2_plus_pes(method_class, **kwargs):
    """
    Calculates the potential energy surface for H2+ for a given method.

    Args:
        method_class: The PySCF calculation class (e.g., scf.UHF, dft.UKS).
        **kwargs: Keyword arguments for the method (e.g., xc='pbe').
    """
    energies = []
    distances = np.arange(0.5, 8.1, 0.25)
    print(f"\n--- Calculating PES for {method_class.__name__} with args: {kwargs} ---")
    print(f"{'Distance (Å)':<15} {'Energy (Hartree)':<20}")
    print("-" * 35)

    for r in distances:
        # Define the molecule: H2+
        # PySCF's gto.Mole() or gto.M() can be used.
        # We set charge=1 and spin=1 (one unpaired electron)
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, r)]
        ]
        mol.basis = 'sto-3g'  # A minimal basis set is sufficient for this demo
        mol.charge = 1
        mol.spin = 1  # 2*S = 1 -> Doublet state
        mol.build()

        # Run the calculation
        # Using .run() is a convenient way to execute the calculation
        if kwargs:
            # For DFT methods that require an 'xc' functional
            mf = method_class(mol, **kwargs).run(conv_tol=1e-7, verbose=0)
        else:
            # For Hartree-Fock
            mf = method_class(mol).run(conv_tol=1e-7, verbose=0)

        energies.append(mf.e_tot)
        print(f"{r:<15.4f} {mf.e_tot:<20.10f}")

    # The exact energy of a separated H atom is -0.5 Hartree
    print("\nNote: The correct dissociation limit (energy at large distance)")
    print("should be the energy of a Hydrogen atom, which is -0.5 Hartree.")
    return distances, energies

# 1. Demonstrate the problem with a DFT functional (PBE)
# UKS = Unrestricted Kohn-Sham (the DFT equivalent of UHF)
# PBE is a common GGA functional that exhibits self-interaction error.
# The energy will incorrectly drop far below -0.5 Hartree.
calculate_h2_plus_pes(dft.UKS, xc='pbe,pbe')

# 2. Demonstrate the fix with Hartree-Fock
# UHF = Unrestricted Hartree-Fock
# HF is free of self-interaction error and will correctly dissociate to -0.5 Hartree.
calculate_h2_plus_pes(scf.UHF)
