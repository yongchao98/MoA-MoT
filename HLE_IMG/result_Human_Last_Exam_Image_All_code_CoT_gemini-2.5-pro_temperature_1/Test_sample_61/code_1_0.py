# Plan:
# 1. Determine the structure of the final product A, which is 2-benzyl-1-indanone.
# 2. Count the number of carbon, hydrogen, and oxygen atoms in this structure.
# 3. Print the resulting molecular formula.

# Atom counts for 2-benzyl-1-indanone:
# Carbons in the indanone core's benzene ring
c_indanone_ring = 6
# Carbons in the indanone core's five-membered ring
c_indanone_five_ring = 3
# Carbons in the benzyl group
c_benzyl = 7
# Total carbon atoms
num_C = c_indanone_ring + c_indanone_five_ring + c_benzyl

# Hydrogens on the indanone core's benzene ring
h_indanone_ring = 4
# Hydrogens on the indanone core's five-membered ring (1 at C-2, 2 at C-3)
h_indanone_five_ring = 3
# Hydrogens on the benzyl group
h_benzyl = 7
# Total hydrogen atoms
num_H = h_indanone_ring + h_indanone_five_ring + h_benzyl

# Oxygen atoms in the ketone group
num_O = 1

# Print the final molecular formula, showing each number as requested.
print(f"The molecular formula of compound A is C{num_C}H{num_H}O{num_O}.")
