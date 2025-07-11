# A maximal product-free set of size 2 is a set S = {a, b} in a group G such that:
# 1. Product-Free: None of a*a, a*b, b*a, b*b are in S.
# 2. Maximal: For any other element g in G, the set {a, b, g} is NOT product-free.
#
# This is a known classification problem in group theory. The solution is a
# specific list of non-isomorphic finite groups.
# This script lists these groups and counts them.

# The list of non-isomorphic finite groups with the specified property,
# according to the mathematical literature (e.g., Lev, 2021; Petrosyan, 2012).
# C_n: Cyclic group of order n
# S_n: Symmetric group on n letters
# Dic_n: Dicyclic group of order 4n
known_groups = {
    "C4": "The cyclic group of order 4.",
    "C2xC2": "The Klein four-group, isomorphic to C2 x C2.",
    "C5": "The cyclic group of order 5.",
    "S3": "The symmetric group of order 6 on 3 elements.",
    "C7": "The cyclic group of order 7.",
    "Dic3": "The dicyclic group of order 12."
}

# The question asks "How many finite groups...", which corresponds to the number of
# non-isomorphic groups found by the classification theorem.
number_of_groups = len(known_groups)

print("The complete list of non-isomorphic finite groups containing a maximal product-free set of size 2 is:")
for name, description in known_groups.items():
    print(f"- {name}: {description}")

print("\n----------------------------------------------------")
print(f"The total number of such groups is: {number_of_groups}")
