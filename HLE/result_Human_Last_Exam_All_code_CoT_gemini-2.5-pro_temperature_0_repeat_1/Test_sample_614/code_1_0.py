# This script solves the riddle by determining the correct order of the numbers 1, 2, 3, 4, and 5.
# The logic is based on the positional clues provided in the riddle.

# "Number 1 is the first..."
# This places the number 1 at the beginning of the sequence.
pos1 = 1

# "...the fifth, who ... lastly follows."
# This places the number 5 at the end of the sequence.
pos5 = 5

# "...three ... protects the last;"
# This means number 3 is in the position just before the last one (the fourth position).
pos4 = 3

# "Number 4 ... always goes behind [two]."
# The remaining numbers are 2 and 4, and the remaining positions are the second and third.
# To satisfy this clue, 2 must come before 4.
pos2 = 2
pos3 = 4

# Now we print the numbers in their correct order to form the final sequence.
# The instruction "output each number in the final equation!" is fulfilled by printing the sequence.
print(f"The final ordered sequence is:")
print(pos1, pos2, pos3, pos4, pos5)