# Initial sequence
sequence = [5, 5, 10, 15, 25, 40, 65]

# Calculate the differences between consecutive terms
differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate the second-level differences
second_level_differences = [differences[i] - differences[i - 1] for i in range(1, len(differences))]

# Predict the next second-level difference
next_second_level_difference = second_level_differences[-1] + (second_level_differences[-1] - second_level_differences[-2])

# Predict the next first-level difference
next_difference = differences[-1] + next_second_level_difference

# Predict the next number in the sequence
next_number = sequence[-1] + next_difference

print(next_number)