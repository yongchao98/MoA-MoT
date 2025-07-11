import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_pes():
    """
    Calculates the potential energy surface of H2+ using
    Restricted Open-shell Hartree-Fock (ROHF) and Unrestricted Hartree-Fock (UHF).
    This demonstrates the failure of the restricted method at dissociation.
    """
    # List of internuclear distances in Angstroms
    distances = np.arange(0.5, 4.1, 0.25)

    print("H2+ Potential Energy Calculation")
    print("Correct dissociation limit should be E(H) + E(H+) = -0.5 Ha")
    print("-" * 60)
    print("{:<15} {:<20} {:<20}".format("Distance (A)", "ROHF Energy (Ha)", "UHF Energy (Ha)"))
    print("-" * 60)

    for r in distances:
        # Define the molecule: H2+
        # Charge = 1, Spin = 1 (one unpaired electron -> doublet)
        mol = gto.Mole()
        mol.atom = [
            ['H', (0., 0., 0.)],
            ['H', (0., 0., r)]
        ]
        mol.charge = 1
        mol.spin = 1
        mol.basis = 'cc-pVDZ'
        mol.build()

        # Restricted Open-shell Hartree-Fock (ROHF) calculation
        mf_rohf = scf.ROHF(mol).run(conv_tol=1e-9, verbose=0)
        e_rohf = mf_rohf.e_tot

        # Unrestricted Hartree-Fock (UHF) calculation
        # This method allows spatial orbitals for alpha and beta electrons to differ,
        # which is how symmetry can be broken.
        mf_uhf = scf.UHF(mol).run(conv_tol=1e-9, verbose=0)
        e_uhf = mf_uhf.e_tot
        
        print("{:<15.4f} {:<20.8f} {:<20.8f}".format(r, e_rohf, e_uhf))

if __name__ == '__main__':
    try:
        calculate_h2_plus_pes()
    except ImportError:
        print("PySCF not found. Please install it to run this script: pip install pyscf")
    except Exception as e:
        print(f"An error occurred: {e}")
