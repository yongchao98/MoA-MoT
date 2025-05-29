# Initial sequence
sequence = [10, 5, 2, 1, 0]

# Calculate differences
differences = [sequence[i] - sequence[i+1] for i in range(len(sequence) - 1)]

# Determine the next difference based on the pattern
next_difference = differences[-1]  # The pattern suggests the next difference is the same as the last one

# Calculate the next term
next_term = sequence[-1] - next_difference

print(next_term)