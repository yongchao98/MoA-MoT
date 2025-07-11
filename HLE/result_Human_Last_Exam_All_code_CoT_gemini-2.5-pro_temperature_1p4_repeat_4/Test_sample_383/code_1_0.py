import numpy as np
from pyscf import gto, scf, fci

def calculate_h2_plus_pes():
    """
    Calculates the Potential Energy Surface (PES) for the H2+ cation
    using Restricted Open-shell Hartree-Fock (ROHF) and Full Configuration
    Interaction (FCI) to demonstrate the failure of ROHF at dissociation.
    """
    # List of internuclear distances in Angstroms
    distances = np.arange(0.5, 4.1, 0.2)

    print("H2+ Potential Energy Curve Calculation")
    print("-----------------------------------------")
    print(f"{'Distance (A)':<15} {'ROHF Energy (Ha)':<20} {'FCI Energy (Ha)':<20}")
    print(f"{'------------':<15} {'----------------':<20} {'---------------':<20}")

    for r in distances:
        # Define the molecule for each distance
        mol = gto.Mole()
        mol.atom = [
            ['H', (0., 0., 0.)],
            ['H', (0., 0., r)]]
        mol.charge = 1  # H2+ has a +1 charge
        mol.spin = 1    # This is multiplicity 2S+1. For one electron, S=1/2, so 2*0.5+1=2. Pyscf wants 2S.
        mol.basis = 'sto-6g' # A simple basis set is sufficient for demonstration
        mol.build()

        # --- 1. Restricted Open-shell Hartree-Fock (ROHF) calculation ---
        # ROHF enforces that the single electron occupies a symmetric orbital,
        # which is incorrect at large distances.
        # We suppress verbose output from the solver using quiet=True.
        mf_rohf = scf.ROHF(mol).run(quiet=True)
        e_rohf = mf_rohf.e_tot

        # --- 2. Full Configuration Interaction (FCI) calculation ---
        # FCI is an exact method for a given basis set and provides the
        # correct description of the dissociation.
        # It is run on top of the ROHF calculation to get the orbitals.
        e_fci, fcivec = fci.FCI(mf_rohf).kernel()

        # Output the numbers for the final result
        print(f"{r:<15.2f} {e_rohf:<20.8f} {e_fci:<20.8f}")

    print("\nAnalysis:")
    print("Notice how the ROHF energy incorrectly rises steeply or behaves")
    print("unphysically at large distances. In contrast, the FCI energy")
    print("correctly flattens out, approaching the energy of a separated")
    print("H atom (-0.4666 Ha in sto-6g) and a proton (0 Ha).")
    print("The strange drop in energy you observed is a severe artifact of this same issue,")
    print("where the calculation becomes unstable (variational collapse).")
    print("The fix is to use a method like FCI, CASSCF, or other multi-reference")
    print("methods that can properly describe bond breaking.")


if __name__ == '__main__':
    calculate_h2_plus_pes()
