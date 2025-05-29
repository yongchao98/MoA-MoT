# Initial sequence
sequence = [-7, -10, -23, -52, -123]

# Calculate first-level differences
first_level_differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Calculate second-level differences
second_level_differences = [first_level_differences[i] - first_level_differences[i-1] for i in range(1, len(first_level_differences))]

# Print the differences
print("First-level differences:", first_level_differences)
print("Second-level differences:", second_level_differences)

# Calculate the next second-level difference assuming a pattern
# Here, we assume the pattern is a linear increase in the second-level differences
# Let's calculate the next second-level difference
next_second_level_difference = second_level_differences[-1] + (second_level_differences[-1] - second_level_differences[-2])

# Calculate the next first-level difference
next_first_level_difference = first_level_differences[-1] + next_second_level_difference

# Calculate the next term in the sequence
next_term = sequence[-1] + next_first_level_difference

# Output the next term
print(next_term)