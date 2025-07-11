# The riddle describes the order of five natural numbers (1, 2, 3, 4, 5).
# Let's determine the order step-by-step based on the clues.

# Clue 1: "the fifth, who never had a thing and lastly follows."
# This means the number 5 is at the end of the sequence.
# Sequence: [_, _, _, _, 5]

# Clue 2: "Number three is the strongest and it protects the last;"
# This means 3 is right before the last number (5).
# Sequence: [_, _, _, 3, 5]

# Clue 3: "Number 1 is the first to have something but is the weakest;"
# This means 1 is at the beginning of the sequence.
# Sequence: [1, _, _, 3, 5]

# Clue 4: "Number 4 likes two the most and it always goes behind."
# This means 4 comes directly after 2, forming a block: (2, 4).
# The only remaining empty spots are the second and third positions.
# The block (2, 4) fits perfectly here.

# Final sequence assembly:
final_sequence = [1, 2, 4, 3, 5]

# The puzzle requests to output each number in the final equation.
# Here is the final ordered sequence.
n1, n2, n3, n4, n5 = final_sequence

print(f"The correct order is: {n1}, {n2}, {n3}, {n4}, {n5}")