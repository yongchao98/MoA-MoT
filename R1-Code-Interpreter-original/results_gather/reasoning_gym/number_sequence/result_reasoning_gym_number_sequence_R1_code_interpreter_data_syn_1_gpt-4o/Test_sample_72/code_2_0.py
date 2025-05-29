# Initial sequence
sequence = [-2, -6, -12, -22, -38, -64, -106, -174]

# Calculate the differences
differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Calculate the second-level differences
second_level_differences = [differences[i] - differences[i-1] for i in range(1, len(differences))]

# Predict the next second-level difference
next_second_level_difference = -42  # Based on the pattern

# Predict the next first-level difference
next_first_level_difference = differences[-1] + next_second_level_difference

# Predict the next term in the sequence
next_term = sequence[-1] + next_first_level_difference

print(next_term)