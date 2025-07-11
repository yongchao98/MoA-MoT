# Based on the riddle's logic, we can deduce the correct order of the five numbers.
# Let's represent the five positions in the sequence as a list.
# sequence = [pos1, pos2, pos3, pos4, pos5]

# 1. "Number 1 is the first to have something but is the weakest;"
# This places the number 1 in the first position.
pos1 = 1

# 2. "...the fifth, who never had a thing and lastly follows."
# This refers to the number 3, placing it in the last (fifth) position.
pos5 = 3

# 3. "Number three is the strongest and it protects the last;"
# The strongest number is 5. "Protects the last" means it's in the position before the last one (the fourth position).
pos4 = 5

# 4. "Number 4 likes two the most and it always goes behind."
# This means the numbers 2 and 4 appear together as a pair: 2, 4.
# The only remaining empty positions are the second and third.
pos2 = 2
pos3 = 4

# Assembling the final sequence
final_sequence = [pos1, pos2, pos3, pos4, pos5]

# Printing the final sequence as requested.
print("The correct order of the numbers is:")
print(*final_sequence, sep=", ")