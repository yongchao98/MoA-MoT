# Initial sequence
sequence = [8, 9, 18, 28, 47]

# Calculate the differences
differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate the second-level differences
second_level_differences = [differences[i] - differences[i - 1] for i in range(1, len(differences))]

# The pattern in second-level differences is alternating between 8 and 1
# Calculate the next difference
next_difference = differences[-1] + 8

# Calculate the next number in the sequence
next_number = sequence[-1] + next_difference

print(next_number)