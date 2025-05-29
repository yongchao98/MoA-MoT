# Initial sequence
sequence = [2, 2, 4, 6, 10, 16]

# Calculate the differences between consecutive terms
differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate the second-level differences
second_level_differences = [differences[i] - differences[i - 1] for i in range(1, len(differences))]

# Predict the next difference
next_difference = differences[-1] + 2  # Assuming the pattern continues

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

print(next_term)