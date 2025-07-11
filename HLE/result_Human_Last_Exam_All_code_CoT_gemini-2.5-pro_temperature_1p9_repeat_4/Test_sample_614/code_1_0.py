# Based on the riddle, the five numbers are determined and placed in order.
# The sequence consists of five natural numbers derived from the clues.

# "Number 1 is the first to have something..." -> The sequence starts with 1.
pos1 = 1

# "...the fifth, who...lastly follows." -> The sequence ends with 5.
pos5 = 5

# "Number three...protects the last;" -> 3 is in the fourth position, before 5.
pos4 = 3

# The remaining numbers are 2 and 4.
# "Number 4...always goes behind [Number 2]." -> 2 must come before 4.
# This places 2 in the second position and 4 in the third.
pos2 = 2
pos3 = 4

# The final ordered sequence is constructed from these positions.
print(f"{pos1}, {pos2}, {pos3}, {pos4}, {pos5}")