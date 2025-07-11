# The reasoning is based on combinatorial analysis rather than direct computation.
# The problem reduces to finding the number of permutations pi' of length m (for m >= 3)
# such that inv(pi')=3 and pi' is 132-avoiding.
# Each such pi' generates a unique permutation of length 333.

# Case m=3:
# Permutations of {1,2,3} with 3 inversions:
# (3,2,1) -> 1 permutation.
# Is (3,2,1) 132-avoiding? Yes.
# Contribution: 1
perm1 = (3, 2, 1)

# Case m=4:
# Permutations of {1,2,3,4} with 3 inversions.
# Let's list them and check for the 132 pattern.
# '2341' -> invs: (2,1),(3,1),(4,1) -> 3. Is it 132-avoiding? Yes.
# '3214' -> invs: (3,2),(3,1),(2,1) -> 3. Is it 132-avoiding? Yes. The disturbed part is {1,2,3}, so this is covered by m=3 case.
# '4123' -> invs: (4,1),(4,2),(4,3) -> 3. Is it 132-avoiding? Yes.
# '1432' -> contains 132 pattern (1,4,3).
# '2413' -> contains 132 pattern (2,4,3).
# '3142' -> contains 132 pattern (1,4,2).
# We are looking for permutations of {1,...,m} that are not permutations of {1,...,m-1}.
# The permutations that are truly of length m=4 are:
perm2 = (2, 3, 4, 1)
perm3 = (4, 1, 2, 3)

# For m >= 5, any permutation with 3 inversions can be shown to contain a 132 pattern.
# Thus, there are no more such permutations.

count = 1 + 1 + 1
print("The number of 1324-avoiding permutations of length 333 with 3 inversions is the sum of contributions from each case.")
print("Case m=3, pi'=(3,2,1): 1")
print("Case m=4, pi'=(2,3,4,1): 1")
print("Case m=4, pi'=(4,1,2,3): 1")
print("Total number = 1 + 1 + 1 = 3")
