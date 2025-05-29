sequence = [8, 5, 3, 2, 2]
differences = [3, 2, 1, 0]

# Calculate the next difference
next_difference = differences[-1] - 1

# Calculate the next term in the sequence
next_term = sequence[-1] + next_difference

print(next_term)