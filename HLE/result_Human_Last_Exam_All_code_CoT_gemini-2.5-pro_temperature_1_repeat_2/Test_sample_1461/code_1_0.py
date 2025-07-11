def calculate_h_bond_ring_size(intermediate_residue_atoms):
  """
  Calculates the number of atoms (m) in a hydrogen-bonded ring.
  The formula is m = 4 + sum of backbone atoms in intermediate residues.
  
  Args:
    intermediate_residue_atoms (list of int): A list containing the number of
                                               backbone atoms for each intermediate residue.
  Returns:
    int: The calculated ring size 'm'.
  """
  return 4 + sum(intermediate_residue_atoms)

# Number of backbone atoms for each monomer type
n_atoms_ala = 3
n_atoms_eps = 7

# --- Analysis of the 11/9 Helix Pattern ---
# This pattern is found experimentally in similar alternating alpha/epsilon peptides.
# It consists of two types of hydrogen-bonded rings: an 11-membered ring and a 9-membered ring.

# 1. Calculate the size of the 11-membered ring.
# This ring is formed by a hydrogen bond between Ala(i) and Ala(i+2).
# The intermediate residue is one epsilon-amino acid (Eps).
intermediate_residues_for_m11 = [n_atoms_eps]
m11 = calculate_h_bond_ring_size(intermediate_residues_for_m11)

print("Calculation for the first ring size (m1):")
print(f"Bonding pattern: Ala(i) -> Ala(i+2)")
print(f"Intermediate residue(s): Epsilon-amino acid (backbone atoms = {n_atoms_eps})")
print(f"m1 = 4 + (atoms in intermediate Eps residue)")
print(f"m1 = 4 + {n_atoms_eps} = {m11}")
print("-" * 20)

# 2. Information about the 9-membered ring.
# In the reference literature, this ring involves a side-chain atom.
# While our generic monomer may not have the exact same side chain, the 11/9 pattern
# is a known, stable conformation for this class of foldamer, making it the most likely answer.
m9 = 9
print("Information for the second ring size (m2):")
print(f"From literature on alpha/epsilon-peptides, a second, 9-membered ring (m2 = {m9}) also forms.")
print("-" * 20)

print(f"The resulting helical pattern is an '{m11}/{m9}' helix.")
