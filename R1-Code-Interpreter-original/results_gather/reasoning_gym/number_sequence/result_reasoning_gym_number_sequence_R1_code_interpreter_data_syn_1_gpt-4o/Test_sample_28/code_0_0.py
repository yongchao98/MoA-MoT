# Initial sequence
sequence = [-6, -14, -28, -50, -86]

# Calculate first-level differences
first_level_differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Calculate second-level differences
second_level_differences = [first_level_differences[i] - first_level_differences[i-1] for i in range(1, len(first_level_differences))]

# Calculate third-level differences
third_level_differences = [second_level_differences[i] - second_level_differences[i-1] for i in range(1, len(second_level_differences))]

# Assume the next third-level difference follows the pattern
next_third_level_difference = third_level_differences[-1] - 4  # Assuming a pattern of decreasing by 4

# Calculate the next second-level difference
next_second_level_difference = second_level_differences[-1] + next_third_level_difference

# Calculate the next first-level difference
next_first_level_difference = first_level_differences[-1] + next_second_level_difference

# Calculate the next term in the sequence
next_term = sequence[-1] + next_first_level_difference

print(next_term)