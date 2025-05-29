# Given sequence
sequence = [4, -1, -2, -8, -15]

# Calculate the differences between consecutive terms
differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate the next difference based on the pattern of differences
# We will assume the pattern in differences continues
next_difference = differences[-1] + (differences[-1] - differences[-2])

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

print(next_term)