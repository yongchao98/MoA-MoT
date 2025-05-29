# Initial sequence
sequence = [8, 9, 18, 28, 47]

# Calculate the differences
differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate the next difference based on the pattern
next_difference = differences[-1] + (differences[-1] - differences[-2])

# Calculate the next number in the sequence
next_number = sequence[-1] + next_difference

print(next_number)