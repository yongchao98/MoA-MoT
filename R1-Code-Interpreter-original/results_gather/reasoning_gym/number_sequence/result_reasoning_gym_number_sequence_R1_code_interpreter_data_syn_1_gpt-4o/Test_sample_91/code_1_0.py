# Initial sequence
sequence = [-9, -14, -21, -32]

# Calculate the differences
differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Calculate the differences of differences
diff_of_diff = [differences[i] - differences[i-1] for i in range(1, len(differences))]

# Calculate the next difference
next_difference = differences[-1] + (diff_of_diff[-1] - 2)

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

print(next_term)