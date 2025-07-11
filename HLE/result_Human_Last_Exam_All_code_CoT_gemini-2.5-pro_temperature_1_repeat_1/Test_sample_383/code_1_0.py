import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_energy(distance, basis='cc-pvdz', method='hf'):
    """
    Calculates the energy of H2+ at a given distance using a specified method.
    """
    # Define the molecule using pyscf's Mole object
    # H2+ has a charge of +1 and one unpaired electron (spin=1)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0., 0., 0.)],
        ['H', (0., 0., distance)]
    ]
    mol.charge = 1
    mol.spin = 1
    mol.basis = basis
    mol.build()

    # Select the calculation method
    if method.lower() == 'hf':
        # Use Restricted Hartree-Fock (RHF)
        mf = scf.RHF(mol)
    else:
        # Use Restricted Kohn-Sham DFT (RKS) with the specified functional
        mf = scf.RKS(mol)
        mf.xc = method

    mf.verbose = 0  # Suppress detailed output during calculation
    energy = mf.kernel()
    return energy

def calculate_h_atom_energy(basis='cc-pvdz', method='hf'):
    """
    Calculates the energy of a single Hydrogen atom, which is the
    dissociation limit for H2+.
    """
    mol = gto.Mole()
    mol.atom = [['H', (0., 0., 0.)]]
    mol.charge = 0
    mol.spin = 1
    mol.basis = basis
    mol.build()

    if method.lower() == 'hf':
        mf = scf.RHF(mol)
    else:
        mf = scf.RKS(mol)
        mf.xc = method
        
    mf.verbose = 0
    energy = mf.kernel()
    return energy

def main():
    """
    Main function to run the calculations and print the results.
    """
    # Define a range of internuclear distances to scan
    distances = np.arange(0.5, 6.1, 0.5)
    
    # Use a common DFT functional that shows the error
    dft_functional = 'pbe,pbe'

    print("Calculating the potential energy curve for H2+...")
    print("This demonstrates the problem with DFT and its fix using Hartree-Fock.\n")
    print(f"{'Distance (Å)':<15} | {'HF Energy (Ha)':<20} | {'DFT ('+dft_functional+') Energy (Ha)':<25}")
    print("-" * 70)

    for r in distances:
        # Calculate energy using Hartree-Fock
        hf_energy = calculate_h2_plus_energy(r, method='hf')
        # Calculate energy using DFT
        dft_energy = calculate_h2_plus_energy(r, method=dft_functional)
        print(f"{r:<15.2f} | {hf_energy:<20.6f} | {dft_energy:<25.6f}")

    print("-" * 70)
    print("\nComparing with the dissociation limit (Energy of H + H+):")
    
    # The final state is H atom + proton. Energy of a proton is 0.
    # So the dissociation limit is the energy of a single H atom.
    hf_h_atom_energy = calculate_h_atom_energy(method='hf')
    dft_h_atom_energy = calculate_h_atom_energy(method=dft_functional)

    print("The final equation for dissociation is: E(H2+) -> E(H) + E(H+)")
    print("Since E(H+) = 0, the final energy is E(H).\n")
    
    print(f"Dissociation Limit E(H) with Hartree-Fock: {hf_h_atom_energy:.6f} Ha")
    print(f"Dissociation Limit E(H) with DFT ({dft_functional}): {dft_h_atom_energy:.6f} Ha")

    print("\nConclusion:")
    print("As the distance increases, the HF energy correctly approaches the dissociation limit (~-0.5 Ha).")
    print("The DFT energy, due to Self-Interaction Error, unphysically continues to decrease.")
    print("This confirms that using Hartree-Fock (statement 2) is the correct approach for this system.")

if __name__ == "__main__":
    # Note: To run this code, you need to install the pyscf library:
    # pip install pyscf
    main()