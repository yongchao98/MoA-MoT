# The set of values for the character corresponding to partition [4,1]
char1_values = [4, 2, 1, 0, 0, -1, -1]

# The set of values for the character corresponding to partition [2,1,1,1]
char2_values = [4, -2, 1, 0, 0, -1, 1]

# Sort each list in ascending order
char1_sorted = sorted(char1_values)
char2_sorted = sorted(char2_values)

# Print the sorted lists, separated by a comma
# This fulfills the request to output each number for each character.
print(f"{char1_sorted},{char2_sorted}")