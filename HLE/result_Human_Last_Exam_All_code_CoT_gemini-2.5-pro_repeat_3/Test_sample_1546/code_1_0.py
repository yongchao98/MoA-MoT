# The irreducible representations of the symmetric group S_n correspond to the partitions of n.
# For S_5, there are 7 partitions. We are interested in the irreducible characters of degree 4.
# These correspond to the partitions [4,1] and its conjugate [2,1,1,1].

# 1. The character chi_[4,1] is the standard character. Its values are calculated as
#    (number of fixed points of a permutation) - 1.
# The 7 conjugacy classes of S_5, represented by their cycle types, are:
# (1,1,1,1,1), (2,1,1,1), (3,1,1), (4,1), (5), (2,2,1), (3,2)
# The number of fixed points for these are, respectively: 5, 3, 2, 1, 0, 1, 0.
# Subtracting 1 gives the character values for chi_[4,1].
char1_values = [4, 2, 1, 0, -1, 0, -1]

# 2. The character chi_[2,1,1,1] is the product of chi_[4,1] and the sign character.
# The signs for the respective conjugacy classes are: +1, -1, +1, -1, +1, +1, -1.
# We multiply the values of char1_values by these signs.
char2_values = [4*1, 2*(-1), 1*1, 0*(-1), (-1)*1, 0*1, (-1)*(-1)]
# This results in the list: [4, -2, 1, 0, -1, 0, 1]

# Sort both lists in ascending order as requested.
char1_sorted = sorted(char1_values)
char2_sorted = sorted(char2_values)

# Print the final lists, separated by a comma.
print(f"{char1_sorted},{char2_sorted}")