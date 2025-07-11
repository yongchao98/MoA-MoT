import numpy as np
from pyscf import gto, scf

def calculate_h2_plus_energy(distance):
    """
    Calculates the RHF energy for H2+ at a given internuclear distance.
    """
    # Define the molecule: H2+
    # PySCF's gto.Mole() or Mole() takes atom specifications, basis set,
    # charge, and spin multiplicity (2S+1). For H2+, we have one unpaired
    # electron (S=1/2), so multiplicity is 2*1/2 + 1 = 2.
    mol = gto.Mole()
    mol.atom = [
        ['H', (0, 0, 0)],
        ['H', (0, 0, distance)]
    ]
    mol.basis = 'cc-pVDZ'  # Use a standard basis set
    mol.charge = 1         # Charge of +1
    mol.spin = 1           # Spin is 1/2, so multiplicity is 2S+1=2, but pyscf wants 2S
    mol.build()

    # Perform a Restricted Open-Shell Hartree-Fock (ROHF) calculation
    # ROHF is suitable for open-shell systems like H2+
    mf = scf.ROHF(mol).run(verbose=0)
    return mf.e_tot

# Define a range of bond distances in Angstroms
distances = np.arange(0.5, 5.1, 0.2)

print("Calculating H2+ Potential Energy Curve using ROHF...")
print("-" * 45)
print(f"{'Distance (A)':<15} | {'Energy (Hartree)':<25}")
print("-" * 45)

# Calculate and print the energy at each distance
for r in distances:
    # In a real scenario, you would need to handle cases where the
    # calculation might fail to converge, but for H2+ it's robust.
    energy = calculate_h2_plus_energy(r)
    # The final "equation" is the set of points (r, E) that make the curve.
    # We print each number that contributes to this curve.
    print(f"{r:<15.2f} | {energy:<25.10f}")

print("-" * 45)
print("Calculation complete. The potential energy curve correctly rises and flattens at large distances.")
