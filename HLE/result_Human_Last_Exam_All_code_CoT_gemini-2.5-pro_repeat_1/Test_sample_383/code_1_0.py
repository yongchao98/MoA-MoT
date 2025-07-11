import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_pes(distances):
    """
    Calculates the potential energy surface for H2+ using Unrestricted Hartree-Fock (UHF).

    UHF is used because it allows for symmetry breaking, which is essential for
    correctly describing the dissociation of H2+ into H + H+. A restricted
    method (ROHF) would incorrectly force the electron to be delocalized over
    both atoms at all distances, leading to an incorrect dissociation energy.

    Args:
        distances (list or np.array): A list of internuclear distances in Angstroms.
    
    Returns:
        A dictionary mapping distance to total energy in Hartrees.
    """
    print("Calculating H2+ Potential Energy Curve using UHF...")
    print("Distance (Angstrom) | Total Energy (Hartree)")
    print("--------------------|-------------------------")
    
    results = {}
    for dist in distances:
        # 1. Define the molecule (H2+)
        # We set charge=1 and spin=1 (one unpaired electron)
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, dist)]
        ]
        mol.charge = 1
        mol.spin = 1
        mol.basis = 'sto-3g' # A minimal basis set for demonstration
        mol.build()

        # 2. Run the Unrestricted Hartree-Fock (UHF) calculation
        # UHF allows the alpha-spin orbital (where the electron is) to localize
        # on one atom at large distances, fixing the dissociation problem.
        mf = scf.UHF(mol)
        
        # Suppress verbose output for a cleaner table
        mf.verbose = 0 
        
        energy = mf.kernel()
        results[dist] = energy
        
        # 3. Print the final result for this distance
        # We output the numbers for the final 'equation' E(R)
        print(f"{dist:<20.4f}| {energy:<25.10f}")
        
    # Find the minimum energy point
    min_dist = min(results, key=results.get)
    min_energy = results[min_dist]
    print("\nCalculation complete.")
    print(f"Equilibrium distance found at approximately {min_dist:.4f} Angstrom")
    print(f"Energy at equilibrium: {min_energy:.10f} Hartree")

    dissoc_dist = max(results.keys())
    dissoc_energy = results[dissoc_dist]
    print(f"Energy at dissociation limit ({dissoc_dist:.1f} Å): {dissoc_energy:.10f} Hartree")

    # The dissociation energy is higher than the equilibrium energy, as expected.
    # The energy of a single H atom in sto-3g is approx -0.46658 Hartree.
    # The dissociation limit of H + H+ should approach this value.
    # A proton (H+) has no electrons and zero energy in this model.
    
if __name__ == '__main__':
    # Define a range of distances to scan, from compressed to dissociated
    bond_distances = np.arange(0.5, 5.5, 0.25)
    calculate_h2_plus_pes(bond_distances)
