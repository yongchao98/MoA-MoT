# Initial sequence
sequence = [-9, -14, -21, -32]

# Calculate the differences between consecutive terms
differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Calculate the second-level differences
second_level_differences = [differences[i] - differences[i-1] for i in range(1, len(differences))]

# Predict the next difference
next_difference = differences[-1] + second_level_differences[-1]

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

print(next_term)