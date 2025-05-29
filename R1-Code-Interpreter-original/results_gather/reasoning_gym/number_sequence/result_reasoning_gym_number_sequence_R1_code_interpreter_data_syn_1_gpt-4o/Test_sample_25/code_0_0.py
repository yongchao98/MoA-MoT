# Initial sequence
sequence = [-1, -4, -8, -15, -26, -44, -73]

# Calculate the differences
differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Calculate the second-level differences
second_level_differences = [differences[i] - differences[i-1] for i in range(1, len(differences))]

# Calculate the next second-level difference
next_second_level_difference = second_level_differences[-1] + second_level_differences[-2]

# Calculate the next first-level difference
next_difference = differences[-1] + next_second_level_difference

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

print(next_term)