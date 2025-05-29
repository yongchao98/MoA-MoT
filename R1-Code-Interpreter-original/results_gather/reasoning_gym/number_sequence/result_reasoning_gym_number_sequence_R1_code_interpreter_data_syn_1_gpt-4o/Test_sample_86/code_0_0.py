# Define the sequence
sequence = [-10, -5, -3, -2, -1, -1]

# Calculate the differences
differences = [sequence[i] - sequence[i-1] for i in range(1, len(sequence))]

# Predict the next difference
next_difference = 0  # Based on the pattern observed

# Calculate the next term
next_term = sequence[-1] + next_difference

# Output the next term
print(next_term)