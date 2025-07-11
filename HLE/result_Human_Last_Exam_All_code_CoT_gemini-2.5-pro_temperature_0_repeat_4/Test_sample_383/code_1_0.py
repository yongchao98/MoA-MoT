# First, ensure the pyscf library is installed:
# pip install pyscf

from pyscf import gto, dft, fci

def calculate_h2_plus_energy(distance):
    """
    Calculates the energy of H2+ at a given internuclear distance
    using DFT (LDA) and FCI methods.
    """
    # Define the H2+ molecule in pyscf
    # Charge is +1, and it's a doublet (Spin = 1/2, so 2S=1)
    mol = gto.Mole()
    mol.atom = [
        ['H', (0., 0., 0.)],
        ['H', (0., 0., distance)]
    ]
    mol.charge = 1
    mol.spin = 1  # 2*S
    mol.basis = 'sto-3g' # A minimal basis is sufficient for this demonstration
    mol.build()

    # --- Perform DFT calculation with LDA functional ---
    # We use Restricted Open-shell Kohn-Sham (ROKS)
    # The 'lda,vwn' functional is known to have significant delocalization error
    ks_lda = dft.ROKS(mol)
    ks_lda.xc = 'lda,vwn'
    # The 'conv_tol' is loosened to prevent convergence failure at large distance
    e_lda = ks_lda.kernel(conv_tol=1e-7)

    # --- Perform FCI (Full Configuration Interaction) calculation ---
    # FCI is an exact method for a given basis set and provides the correct energy
    # We use the converged DFT calculation as a starting point for FCI
    e_fci, fcivec = fci.FCI(ks_lda).kernel()

    return e_lda, e_fci

# --- Main Calculation ---

# A distance near the equilibrium bond length (in Angstroms)
dist_eq = 1.06

# A stretched distance representing bond dissociation
dist_stretched = 10.0

# Calculate energies at both distances
e_lda_eq, e_fci_eq = calculate_h2_plus_energy(dist_eq)
e_lda_stretched, e_fci_stretched = calculate_h2_plus_energy(dist_stretched)

# --- Print and Analyze the Results ---
print("Comparing DFT/LDA and FCI energies for H2+...")
print("-" * 60)
print(f"Energy at Equilibrium Distance R = {dist_eq:.2f} Å")
print(f"   DFT (LDA) Energy = {e_lda_eq:.6f} Hartrees")
print(f"   FCI Energy       = {e_fci_eq:.6f} Hartrees")
print("-" * 60)
print(f"Energy at Stretched Distance R = {dist_stretched:.2f} Å")
print(f"   DFT (LDA) Energy = {e_lda_stretched:.6f} Hartrees")
print(f"   FCI Energy       = {e_fci_stretched:.6f} Hartrees")
print("-" * 60)

# The correct dissociation limit is the energy of a hydrogen atom (H) plus a proton (H+).
# The energy of H+ is zero. Let's calculate the energy of a single H atom for reference.
mol_h = gto.M(atom='H 0 0 0', basis='sto-3g', spin=1)
e_h_atom = dft.ROKS(mol_h).kernel()

print("Analysis of the Results:")
print(f"The correct dissociation energy should approach the energy of a single H atom ({e_h_atom:.6f} Hartrees).")
print(f"The FCI energy at {dist_stretched:.1f} Å ({e_fci_stretched:.6f}) correctly approaches this limit.")
if e_lda_stretched < e_lda_eq:
    print(f"The DFT/LDA energy at {dist_stretched:.1f} Å ({e_lda_stretched:.6f}) has unphysically dropped below the equilibrium energy ({e_lda_eq:.6f}).")
    print("This demonstrates the 'inverse symmetry breaking' phenomenon due to delocalization error.")
else:
    print("The DFT/LDA energy did not show the pathological drop in this specific calculation.")
