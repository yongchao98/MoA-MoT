# Original sequence
sequence = [-5, -5, -10, -15, -25, -40]

# Calculate the differences
differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Calculate the next difference
next_difference = differences[-1] - 5

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

print(next_term)