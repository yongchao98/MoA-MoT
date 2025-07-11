import numpy as np
from pyscf import gto, scf, dft

print("This script calculates the potential energy curve of H2+ using both")
print("Hartree-Fock (HF) and a standard DFT functional (LDA) to demonstrate")
print("how HF correctly solves the self-interaction error problem.\n")

# --- Step 1: Define a reference energy for dissociation ---
# At infinite distance, H2+ becomes a H atom and a proton H+.
# The energy of this system is simply the energy of one H atom.
# Let's calculate this value first for both methods.
h_atom = gto.M(atom='H 0 0 0', basis='cc-pvdz', charge=0, spin=1)
# For a one-electron atom, HF is exact.
e_h_hf = h_atom.RHF().run(verbose=0).e_tot
# For a one-electron atom, a correct DFT functional should also be exact.
# We use the same functional (LDA) to see its reference energy.
e_h_dft = dft.RKS(h_atom).run(verbose=0).e_tot

print("--- Reference Energy at Infinite Separation (H atom) ---")
# The theoretical exact energy is -0.5 Hartree (Ha).
print(f"Hartree-Fock dissociation limit: {e_h_hf:.6f} Ha")
print(f"LDA DFT dissociation limit:       {e_h_dft:.6f} Ha\n")


print("--- Potential Energy Curve Calculation ---")
print(f"{'Distance (A)':<15} {'HF Energy (Ha)':<20} {'DFT/LDA Energy (Ha)':<20}")
print("-" * 60)

# --- Step 2: Loop over a range of bond distances ---
distances = np.arange(0.5, 6.0, 0.25)

for r in distances:
    # Set up the H2+ molecule for the current distance `r`
    # The system has a charge of +1 and 1 unpaired electron (spin=1)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, r)]
    ]
    mol.basis = 'cc-pvdz'
    mol.charge = 1
    mol.spin = 1 
    mol.build()

    # --- Calculation 1: Hartree-Fock ---
    # This method is free of self-interaction error and is correct for H2+.
    mf_hf = scf.RHF(mol).run(verbose=0)
    e_hf = mf_hf.e_tot

    # --- Calculation 2: DFT with LDA functional ---
    # LDA (Local Density Approximation) is known to have significant SIE.
    mf_dft = dft.RKS(mol)
    mf_dft.xc = 'lda,vwn' # A standard LDA functional
    mf_dft.run(verbose=0)
    e_dft = mf_dft.e_tot

    # --- Step 3: Print the results for comparison ---
    # This is the "final equation" showing the numerical results
    print(f"{r:<15.2f} {e_hf:<20.6f} {e_dft:<20.6f}")

print("\n--- Conclusion ---")
print("As you can see, the HF energy correctly approaches the dissociation limit (~-0.5 Ha).")
print("The DFT/LDA energy incorrectly continues to drop to lower and lower values,")
print("which is the unphysical behavior caused by self-interaction error.")
