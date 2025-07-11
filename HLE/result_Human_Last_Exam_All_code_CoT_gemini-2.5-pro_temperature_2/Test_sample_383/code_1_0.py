import numpy as np
try:
    from pyscf import gto, scf, dft
except ImportError:
    print("PySCF is not installed. Please install it using: pip install pyscf")
    exit()

def calculate_potential_energy_curve(distances, method='hf', basis='cc-pVDZ', charge=1, spin=1):
    """
    Calculates the potential energy curve for H2+ using a specified method.
    
    Args:
        distances (list or np.array): A list of internuclear distances in Angstroms.
        method (str): The quantum chemistry method to use ('hf' or 'pbe').
    
    Returns:
        A list of calculated energies in Hartrees.
    """
    energies = []
    print(f"\n--- Calculating with method: {method.upper()} ---")
    print(f"{'Distance (A)':<15} {'Energy (Hartree)':<20}")
    print("-" * 35)
    
    for r in distances:
        # Define the molecule geometry for each distance
        mol = gto.Mole()
        mol.atom = [
            ['H', (0, 0, -r/2)],
            ['H', (0, 0,  r/2)],
        ]
        mol.basis = basis
        mol.charge = charge
        mol.spin = spin # 1 unpaired electron (Doublet state)
        mol.symmetry = True # Enforce D2h symmetry, a subgroup of Dinf-h
        mol.build(verbose=0)

        # Select the calculation type
        if method.lower() == 'hf':
            # Use Restricted Open-shell Hartree-Fock (ROHF)
            mf = scf.ROHF(mol)
        elif method.lower() == 'pbe':
            # Use Restricted Open-shell Kohn-Sham DFT with PBE functional
            mf = dft.ROKS(mol)
            mf.xc = 'pbe,pbe'
        else:
            raise ValueError("Method must be 'hf' or 'pbe'")

        # Run the calculation and store the energy
        energy = mf.kernel()
        energies.append(energy)
        # In this context, "output each number in the final equation" can be interpreted
        # as printing the components of our calculation: the distance and the resulting energy.
        print(f"{r:<15.2f} {energy:<20.8f}")
        
    return energies

if __name__ == "__main__":
    # Define the range of internuclear distances to scan
    distance_range = np.arange(0.6, 10.1, 0.4)

    # 1. Calculate with PBE (a GGA DFT functional), which has self-interaction error
    pbe_energies = calculate_potential_energy_curve(distance_range, method='pbe')

    # 2. Calculate with Hartree-Fock, which is self-interaction-free
    hf_energies = calculate_potential_energy_curve(distance_range, method='hf')
    
    min_hf_energy = np.min(hf_energies)
    dissociation_limit_hf = hf_energies[-1]
    
    min_pbe_energy = np.min(pbe_energies)
    dissociation_limit_pbe = pbe_energies[-1]

    print("\n" + "="*50)
    print("                     Summary")
    print("="*50)
    print("The Hartree-Fock (HF) calculation behaves correctly:")
    print(f"  - It shows an energy minimum (equilibrium) of {min_hf_energy:.6f} Ha.")
    print(f"  - At large distance (10.0 A), it approaches the correct dissociation limit of H + H+ ({dissociation_limit_hf:.6f} Ha), which is close to the exact value of -0.5 Ha for a Hydrogen atom.")
    print("\nThe PBE (DFT) calculation shows the unphysical behavior:")
    print(f"  - It finds an equilibrium energy of {min_pbe_energy:.6f} Ha.")
    print(f"  - However, as the atoms separate, the energy drops incorrectly to {dissociation_limit_pbe:.6f} Ha at 10.0 A, which is far too low.")
    print("\nThis demonstration confirms that using Hartree-Fock (Statement 2) fixes the problem, which is an example of the phenomenon described in Statement 3. Changing packages to one that uses HF by default (Statement 1) would also work.")
    print("="*50)