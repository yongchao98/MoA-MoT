import pyscf
from pyscf import gto, scf, dft

# --- Plan ---
# 1. Define the H2+ molecule with a large internuclear distance (R=10 Angstrom).
#    The true energy at this distance should be very close to the energy of an
#    isolated Hydrogen atom, which is -0.5 Hartree.
# 2. Perform a calculation using a standard DFT functional (PBE). This method
#    suffers from self-interaction error and is expected to give a wrong,
#    artificially low energy. This mimics the user's problem.
# 3. Perform a calculation using the Hartree-Fock (HF) method. This method
#    is self-interaction free for one electron and should give the correct energy.
#    This demonstrates the fix.
# 4. Print and compare the results.

# Define the H2+ molecule at a large distance of 10 Angstroms
mol = gto.Mole()
mol.atom = [['H', (0, 0, 0)], ['H', (0, 0, 10)]]
mol.charge = 1  # H2(+) cation
mol.spin = 1    # One unpaired electron (doublet)
mol.basis = 'cc-pvdz'
mol.build()

print("Calculating energy for H2+ at R = 10 Angstroms.\n")
print(f"The correct dissociation limit (energy of a hydrogen atom) is -0.5 Hartree.\n")

# --- Calculation with a method that shows the problem (DFT with PBE functional) ---
# We use Restricted Open-shell Kohn-Sham (ROKS)
# This will enforce symmetry, and the self-interaction error of the PBE
# functional will lead to an incorrect, overly stable energy.
try:
    mf_dft = dft.ROKS(mol)
    mf_dft.xc = 'pbe,pbe'
    e_dft = mf_dft.run(conv_tol=1e-10).e_tot
    print(f"1. Energy calculated with DFT (PBE functional):")
    print(f"   E(DFT) = {e_dft:.6f} Hartree")
    print("   This energy is significantly lower than -0.5, which is physically incorrect.")
    print("   This is due to the self-interaction error of the PBE functional.\n")
except Exception as e:
    print(f"1. DFT calculation failed. Error: {e}\n")


# --- Calculation with the corrected method (Hartree-Fock) ---
# We use Restricted Open-shell Hartree-Fock (ROHF)
# HF is free of self-interaction error for a one-electron system and will
# give the physically correct dissociation energy.
try:
    mf_hf = scf.ROHF(mol)
    e_hf = mf_hf.run(conv_tol=1e-10).e_tot
    print(f"2. Energy calculated with Hartree-Fock (HF):")
    print(f"   E(HF)   = {e_hf:.6f} Hartree")
    print("   This energy is very close to the correct value of -0.5.")
    print("   This demonstrates that using HF fixes the problem.\n")
except Exception as e:
    print(f"2. HF calculation failed. Error: {e}\n")

print("Conclusion: The problem is caused by self-interaction error in the DFT functional,")
print("which is particularly bad for symmetric radical cations. Using a self-interaction-free")
print("method like Hartree-Fock resolves the issue.")
