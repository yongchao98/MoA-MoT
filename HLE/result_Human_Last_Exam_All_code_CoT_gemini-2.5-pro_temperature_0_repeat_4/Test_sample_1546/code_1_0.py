# The irreducible characters of degree 4 for S5 correspond to the partitions [4,1] and [2,1,1,1].
# The values for the 7 conjugacy classes are calculated as described in the plan.

# Character values for the partition [4,1]
char1_values = [4, 2, 0, 1, -1, 0, -1]

# Character values for the partition [2,1,1,1]
char2_values = [4, -2, 0, 1, 1, 0, -1]

# Sort each list of character values in ascending order
sorted_char1 = sorted(char1_values)
sorted_char2 = sorted(char2_values)

# Print the final sorted lists, separated by a comma
print(f"{sorted_char1},{sorted_char2}")