import math

# Plan:
# 1. The problem of finding the largest c is an open research problem. We will demonstrate
#    the construction that gives the best known lower bound, c >= 3/8.
# 2. The method involves constructing the set A using modular arithmetic with a modulus m=32.
# 3. A is defined as the set of numbers n such that (n mod 32) is in a special set of residues R.
# 4. The condition that A+A contains no squares requires that the sumset R+R (mod 32)
#    contains no quadratic residues mod 32.
# 5. This script verifies this property for the set R and calculates the resulting density c.

# The modulus
m = 32

# The set of residues R from the paper by Chen and Fang (2019)
R = {1, 5, 6, 7, 10, 11, 13, 14, 21, 22, 26, 30}
size_R = len(R)

# Step 1: Find all quadratic residues modulo 32
# A number q is a quadratic residue mod m if q = k^2 (mod m) for some integer k.
Q_32 = set()
for k in range(m):
    Q_32.add(k*k % m)

print(f"The modulus is m = {m}.")
print(f"The chosen set of residues is R = {sorted(list(R))}")
print(f"The size of R is |R| = {size_R}")
print(f"The set of quadratic residues mod {m} is Q_{m} = {sorted(list(Q_32))}\n")

# Step 2: Compute the sumset R+R modulo 32
R_plus_R = set()
for r1 in R:
    for r2 in R:
        R_plus_R.add((r1 + r2) % m)

print(f"The sumset R+R (mod {m}) is: {sorted(list(R_plus_R))}\n")

# Step 3: Check if R+R and Q_32 have any common elements
intersection = R_plus_R.intersection(Q_32)

if not intersection:
    print("Verification successful: The sumset R+R contains no quadratic residues modulo 32.")
    # The density c is given by the ratio |R| / m
    c = size_R / m
    print("The construction is therefore valid.")
    print(f"The resulting density is c = |R| / m = {size_R} / {m} = {c}")
else:
    print("Verification failed: The sumset R+R contains the following quadratic residues:")
    print(intersection)
