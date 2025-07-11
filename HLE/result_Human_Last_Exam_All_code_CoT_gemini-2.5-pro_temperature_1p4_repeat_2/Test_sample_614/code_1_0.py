# The riddle describes a specific order for the numbers 1, 2, 3, 4, and 5.
# Based on the clues, the sequence is determined as follows:
# - The fifth (5) "lastly follows", so it is in the last position.
# - Number three (3) "protects the last", so it is in the fourth position.
# - Number one (1) is "the first", so it is in the first position.
# - Number four (4) "goes behind" number two (2), filling the remaining second and third spots.

# The final sequence of numbers.
final_order = [1, 2, 4, 3, 5]

# The prompt asks to "output each number in the final equation!".
# We will print each number in the solved sequence.
print("The correct order is:")
print(final_order[0], final_order[1], final_order[2], final_order[3], final_order[4])
