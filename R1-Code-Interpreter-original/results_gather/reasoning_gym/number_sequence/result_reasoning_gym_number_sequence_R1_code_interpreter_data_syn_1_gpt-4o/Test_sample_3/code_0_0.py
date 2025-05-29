# Given sequence
sequence = [-7, -1, -5, -1, -3, 0, 0, 3]

# Calculate differences
differences = [sequence[i+1] - sequence[i] for i in range(len(sequence) - 1)]

# Predict the next difference
# Based on the pattern, we assume the next difference might be -2
next_difference = -2

# Calculate the next term
next_term = sequence[-1] + next_difference

# Output the next term
print(next_term)