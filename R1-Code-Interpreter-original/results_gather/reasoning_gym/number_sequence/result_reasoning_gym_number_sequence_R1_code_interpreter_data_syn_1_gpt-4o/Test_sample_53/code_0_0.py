# Initial sequence
sequence = [-3, -3, -6, -9, -15]

# Calculate the differences
differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Calculate the next difference based on the pattern
next_difference = -3  # Following the pattern of alternating -3 and 0

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

# Output the next term
print(next_term)