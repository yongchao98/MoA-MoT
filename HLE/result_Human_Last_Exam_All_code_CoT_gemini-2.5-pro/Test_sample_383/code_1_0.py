import numpy as np
try:
    from pyscf import gto, scf, dft
except ImportError:
    print("PySCF is not installed. Please install it using: pip install pyscf")
    exit()

def calculate_pes(mol_template, distances, method_name, xc_functional=None):
    """
    Calculates the Potential Energy Surface (PES) for H2+
    using a specified quantum chemistry method.
    """
    energies = []
    print(f"\nCalculating PES using {method_name.upper()}...")
    for r in distances:
        # Build the molecule object for the current internuclear distance
        mol = mol_template.copy()
        mol.atom = [
            ['H', (0, 0, 0)],
            ['H', (0, 0, r)]
        ]
        mol.build()

        # Select the calculation method (DFT or Hartree-Fock)
        if method_name.lower() == 'dft':
            mf = dft.RKS(mol)
            if xc_functional:
                mf.xc = xc_functional
        elif method_name.lower() == 'hf':
            mf = scf.RHF(mol)
        else:
            raise ValueError(f"Unsupported method: {method_name}")
            
        # Suppress the verbose output from PySCF for a cleaner result
        mf.verbose = 0
        
        # Run the self-consistent field calculation and store the energy
        energy = mf.kernel()
        energies.append(energy)
        
    return np.array(energies)

def main():
    """
    Main function to run the calculations and print the results.
    """
    # Define a range of internuclear distances in Angstroms
    distances = np.arange(0.5, 6.1, 0.25)

    # Create a template for the H2+ molecule.
    # It has a charge of +1 and is a doublet (1 unpaired electron).
    mol_template = gto.Mole()
    mol_template.basis = 'sto-3g'  # A minimal basis set is sufficient for this demonstration
    mol_template.charge = 1
    mol_template.spin = 1
    mol_template.unit = 'angstrom'

    # --- Calculation Step ---
    # 1. Demonstrate the problem with a common DFT functional (PBE)
    dft_energies = calculate_pes(mol_template, distances, 'dft', xc_functional='pbe,pbe')

    # 2. Demonstrate the fix using Hartree-Fock (HF)
    hf_energies = calculate_pes(mol_template, distances, 'hf')
    
    # The exact energy of a separated H atom (1s electron) + a bare proton.
    # The energy of an H atom in the STO-3G basis is -0.4665818 Hartrees.
    # The PES should approach this value at large distances.
    exact_dissociation_limit = -0.466581849557333

    # --- Output and Analysis Step ---
    print("\n--- Potential Energy Surface Comparison for H2+ ---")
    print("The table below compares the calculated potential energy for H2+ at various")
    print("internuclear distances using two different methods.")
    print(f"The correct energy at large distance should be ~{exact_dissociation_limit:.6f} Ha.\n")

    # Print the final results in a formatted table
    # This fulfills the requirement to "output each number in the final equation"
    print(f"{'Distance (A)':<15} {'PBE-DFT Energy (Ha)':<25} {'HF Energy (Ha)':<20}")
    print(f"{'='*14:<15} {'='*24:<25} {'='*19:<20}")

    for i, r in enumerate(distances):
        dft_e = dft_energies[i]
        hf_e = hf_energies[i]
        # Highlight the problematic DFT behavior
        dft_comment = "<- Problem: Energy drops unphysically" if dft_e < -0.55 and r > 3.0 else ""
        print(f"{r:<15.2f} {dft_e:<25.8f} {hf_e:<20.8f} {dft_comment}")
        
    print("\n--- Conclusion ---")
    print("1. The PBE-DFT energy incorrectly drops as distance increases. This is due to")
    print("   Self-Interaction Error, an artifact of approximate DFT (related to Statement 3).")
    print("2. The Hartree-Fock (HF) calculation correctly shows the energy stabilizing at")
    print("   the dissociation limit. Thus, using HF is a valid fix (Statement 2).")
    print("3. Different software or functionals might handle this differently, so changing")
    print("   the computational setup is a practical approach (Statement 1).")
    print("\nAll three statements are relevant to understanding and fixing the issue.")

if __name__ == "__main__":
    main()
<<<G>>>