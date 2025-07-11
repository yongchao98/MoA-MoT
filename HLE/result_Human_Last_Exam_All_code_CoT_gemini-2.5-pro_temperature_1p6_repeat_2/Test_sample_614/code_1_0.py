# Let's determine the order of the five numbers based on the riddle's clues.
# The numbers to arrange are 1, 2, 3, 4, 5.

# Clue: "Number 1 is the first to have something..."
# This places 1 in the first position.
pos1 = 1

# Clue: "...the fifth, who never had a thing and lastly follows."
# This places 5 in the last (5th) position.
pos5 = 5

# Clue: "Number three is the strongest and it protects the last;"
# This means 3 is in the position right before the last one (the 4th position).
pos4 = 3

# Clue: "Number 4 likes two the most and it always goes behind."
# This means the order of 2 and 4 must be '2, 4'.
# The only remaining positions are the 2nd and 3rd.
# Therefore, 2 must be in the 2nd position and 4 in the 3rd.
pos2 = 2
pos3 = 4

# The final sequence is constructed from these deductions.
final_sequence = [pos1, pos2, pos3, pos4, pos5]

# Print the final sequence.
print("The correct order of the numbers is:")
print(f"{final_sequence[0]}, {final_sequence[1]}, {final_sequence[2]}, {final_sequence[3]}, {final_sequence[4]}")