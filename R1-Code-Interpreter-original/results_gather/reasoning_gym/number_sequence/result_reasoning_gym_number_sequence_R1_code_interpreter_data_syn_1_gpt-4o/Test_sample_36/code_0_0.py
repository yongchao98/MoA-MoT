# Define the sequence
sequence = [-5, -3, -2, -1, -1, -1, -1]

# Check the differences between consecutive numbers
differences = [sequence[i+1] - sequence[i] for i in range(len(sequence) - 1)]

# Print the differences
print("Differences:", differences)

# Predict the next number based on the pattern
next_number = sequence[-1]  # Since the sequence stabilizes at -1
print("Next number in the sequence:", next_number)