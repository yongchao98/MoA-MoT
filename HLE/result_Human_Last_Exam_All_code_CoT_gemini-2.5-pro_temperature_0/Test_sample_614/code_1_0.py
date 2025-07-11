# Based on the riddle's clues, the sequence is determined as follows:
# 1. "Number 1 is the first..." -> Position 1 is 1.
# 2. "...the fifth, who...lastly follows." -> Position 5 is 5.
# 3. "Number three...protects the last;" -> Position 4 is 3.
# 4. "Number 4...always goes behind [Number 2]." -> The pair 2, 4 fills the remaining spots.
# This results in the sequence: 1, 2, 4, 3, 5.

final_sequence = [1, 2, 4, 3, 5]

# Print each number in the final sequence
print("The correct order of the numbers is:")
print(*final_sequence, sep=", ")