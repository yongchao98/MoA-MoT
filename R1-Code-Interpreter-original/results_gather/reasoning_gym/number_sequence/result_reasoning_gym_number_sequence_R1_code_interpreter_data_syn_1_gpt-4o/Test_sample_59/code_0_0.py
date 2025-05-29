# Initial sequence
sequence = [-4, -4, -8, -12, -20, -32, -52]

# Calculate first-level differences
first_level_differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate second-level differences
second_level_differences = [first_level_differences[i] - first_level_differences[i - 1] for i in range(1, len(first_level_differences))]

# Predict the next second-level difference
next_second_level_difference = second_level_differences[-1] - 4

# Predict the next first-level difference
next_first_level_difference = first_level_differences[-1] + next_second_level_difference

# Predict the next term in the sequence
next_term = sequence[-1] + next_first_level_difference

print(next_term)