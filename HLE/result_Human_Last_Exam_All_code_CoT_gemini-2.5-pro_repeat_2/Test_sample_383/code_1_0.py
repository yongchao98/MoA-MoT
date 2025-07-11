import numpy as np
from pyscf import gto, scf, dft

def calculate_h2_plus_energy(distance):
    """
    Calculates the energy of H2+ at a given internuclear distance
    using both Hartree-Fock and DFT (PBE).
    """
    # Define the molecule
    mol = gto.Mole()
    # Set up the geometry, charge, and spin
    mol.atom = [
        ['H', (0., 0., 0.)],
        ['H', (0., 0., distance)]]
    mol.charge = 1  # H2(+)
    mol.spin = 1    # One unpaired electron (doublet)
    mol.basis = 'cc-pvdz'
    mol.build()

    # --- Hartree-Fock Calculation (Correct for H2+) ---
    # For a one-electron system, RHF/ROHF is exact and free of self-interaction error.
    mf_hf = scf.RHF(mol)
    # Suppress verbose output for cleaner printing
    mf_hf.verbose = 0
    e_hf = mf_hf.kernel()

    # --- DFT Calculation with PBE functional (Shows the error) ---
    # PBE is a common GGA functional that suffers from self-interaction error.
    mf_dft = dft.RKS(mol)
    mf_dft.xc = 'pbe,pbe'
    # Suppress verbose output
    mf_dft.verbose = 0
    e_dft = mf_dft.kernel()

    return e_hf, e_dft

def main():
    """
    Main function to calculate and print the potential energy surface.
    """
    print("Calculating the H2+ potential energy curve with HF and DFT (PBE)...")
    print("Note the incorrect dip in the PBE energy at large distances.")
    print("-" * 65)
    print(f"{'Distance (A)':<15} {'HF Energy (Ha)':<25} {'PBE Energy (Ha)':<25}")
    print("-" * 65)

    # Define the range of distances to scan
    distances = np.arange(0.5, 8.1, 0.5)
    
    # Equilibrium distance for H2+ is around 1.06 A
    equilibrium_dist = 1.06
    
    # Calculate and print energy at equilibrium distance first for reference
    e_hf_eq, e_dft_eq = calculate_h2_plus_energy(equilibrium_dist)
    print(f"{equilibrium_dist:<15.4f} {e_hf_eq:<25.10f} {e_dft_eq:<25.10f} <-- Equilibrium")

    # Calculate for the rest of the distances
    for r in distances:
        if np.isclose(r, equilibrium_dist): # Avoid recalculating
            continue
        e_hf, e_dft = calculate_h2_plus_energy(r)
        print(f"{r:<15.4f} {e_hf:<25.10f} {e_dft:<25.10f}")

    print("-" * 65)
    print("As distance -> infinity, the correct energy should be -0.5 Ha (isolated H atom).")
    print("The HF energy correctly approaches this limit from below.")
    print("The PBE energy incorrectly drops far below the equilibrium energy.")


if __name__ == '__main__':
    main()
