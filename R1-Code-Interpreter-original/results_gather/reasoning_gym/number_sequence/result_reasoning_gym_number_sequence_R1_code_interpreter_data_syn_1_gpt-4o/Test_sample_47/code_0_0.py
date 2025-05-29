# Define the sequence
sequence = [-8, -8, -16, -24, -40, -64, -104, -168]

# Calculate the first-level differences
first_level_differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate the second-level differences
second_level_differences = [first_level_differences[i] - first_level_differences[i - 1] for i in range(1, len(first_level_differences))]

# Calculate the next second-level difference
next_second_level_difference = second_level_differences[-1] - 8

# Calculate the next first-level difference
next_first_level_difference = first_level_differences[-1] + next_second_level_difference

# Calculate the next term in the sequence
next_term = sequence[-1] + next_first_level_difference

# Output the next term
print(next_term)