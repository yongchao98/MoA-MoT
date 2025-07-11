# This script solves the riddle by deducing the order of the numbers based on the provided clues.

# We need to order five numbers. Let's start with an empty list of size 5.
# We will fill it according to the clues.
sequence = [None] * 5

# Clue: "...the fifth, who never had a thing and lastly follows."
# This identifies the fifth number as 0 and places it at the last position (index 4).
sequence[4] = 0

# Clue: "Number three is the strongest and it protects the last;"
# "Protects the last" means it comes just before the last number. So, 3 is at index 3.
sequence[3] = 3

# Clue: "Number 1 is the first to have something..."
# This places 1 at the first position (index 0).
sequence[0] = 1

# Clue: "Number 4 likes two the most and it always goes behind."
# This means 4 always follows 2, forming a pair (2, 4).
# The remaining empty slots are at index 1 and 2. The pair (2, 4) fits perfectly.
sequence[1] = 2
sequence[2] = 4

# Now the sequence is complete. Let's print the result.
# The riddle asks to "Put the numbers in the correct order", which we represent as a final sequence.
print("The final sequence of numbers is:")
print(f"{sequence[0]}, {sequence[1]}, {sequence[2]}, {sequence[3]}, {sequence[4]}")