# The irreducible characters of degree 4 for the symmetric group S_5 correspond to the partitions [4,1] and [2,1,1,1].
# The values of these characters on the 7 conjugacy classes of S_5 have been pre-calculated.
# This script will store those values in lists, sort them, and print them in the requested format.

# Character values corresponding to partition [4,1]
char1_values = [4, 2, 1, 0, -1, 0, -1]
# Character values corresponding to partition [2,1,1,1]
char2_values = [4, -2, 1, 0, -1, 0, 1]

# Sort the character values in ascending order
char1_values.sort()
char2_values.sort()

# Print the two sorted lists, separated by a comma, as specified.
# The str() function is used to get the exact string representation of the list.
print(str(char1_values) + ", " + str(char2_values))