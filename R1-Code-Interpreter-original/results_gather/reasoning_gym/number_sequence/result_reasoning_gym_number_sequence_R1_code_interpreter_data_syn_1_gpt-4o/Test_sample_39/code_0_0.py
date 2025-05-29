# Initial sequence
sequence = [-6, -6, -12, -18, -30, -48, -78]

# Calculate the differences
differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate the next difference
next_difference = differences[-1] + differences[-2]

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

# Output the next term
print(next_term)