import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_pes(distances):
    """
    Calculates and prints the potential energy curve for the H2+ molecular cation
    using both a DFT functional (PBE) known to have self-interaction error and
    Hartree-Fock (HF), which is self-interaction free for one-electron systems.

    This demonstrates the failure of the DFT functional at long bond distances
    and the correct behavior of the HF method.
    """
    print("Calculating Potential Energy Curve for H2+")
    print("The DFT (PBE) energy will show an unphysical drop at large distances,")
    print("while the Hartree-Fock (HF) energy will correctly approach the dissociation limit.")
    print("-" * 65)
    print("{:<15s} {:<25s} {:<25s}".format("Distance (A)", "DFT Energy (PBE) [Hartree]", "HF Energy [Hartree]"))
    print("-" * 65)

    # Energy of a single H atom with this basis is the correct dissociation limit for H2+
    h_atom_mol = gto.M(atom='H', basis='sto-3g', charge=0, spin=1)
    dissociation_limit = scf.UHF(h_atom_mol).run().e_tot

    for r in distances:
        # Define the H2+ molecule for the current distance
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, r)]
        ]
        mol.basis = 'sto-3g'  # Minimal basis set is sufficient for demonstration
        mol.charge = 1        # Cation
        mol.spin = 1          # One unpaired electron (Doublet state)
        mol.build(verbose=0)  # Suppress verbose build output

        # --- DFT Calculation (using PBE, a GGA functional) ---
        mf_pbe = scf.UKS(mol)
        mf_pbe.xc = 'pbe,pbe' # Use the PBE functional
        try:
            # We add a check for convergence failure at long distances
            mf_pbe.kernel()
            if mf_pbe.converged:
                e_pbe = mf_pbe.e_tot
            else:
                e_pbe = float('nan')
        except Exception:
            e_pbe = float('nan') # Mark as calculation failed

        # --- Hartree-Fock Calculation ---
        mf_hf = scf.UHF(mol)
        try:
            mf_hf.kernel()
            if mf_hf.converged:
                e_hf = mf_hf.e_tot
            else:
                e_hf = float('nan')
        except Exception:
            e_hf = float('nan')

        # Print the results for the current distance in a formatted table
        print("{:<15.2f} {:<25.8f} {:<25.8f}".format(r, e_pbe, e_hf))

    print("-" * 65)
    print(f"Correct Dissociation Limit (Energy of H atom): {dissociation_limit:.8f} Hartree")

if __name__ == '__main__':
    # Define a range of bond distances in Angstroms to scan
    # We start near the equilibrium distance (~1.06 A) and stretch the bond
    distances_to_scan = np.arange(0.8, 10.1, 0.4)
    calculate_h2_plus_pes(distances_to_scan)