# Given sequence
sequence = [3, 2, 4, 5]

# Calculate differences between consecutive numbers
differences = [sequence[i+1] - sequence[i] for i in range(len(sequence) - 1)]

# Print the differences
print("Differences:", differences)

# Predict the next number in the sequence
# Assuming a pattern in differences, let's try to find the next difference
# and compute the next number in the sequence

# Let's assume the differences repeat or follow a pattern
# We can try a simple approach by repeating the differences
next_difference = differences[len(sequence) % len(differences)]

# Calculate the next number
next_number = sequence[-1] + next_difference

# Output the next number
print("Next number in the sequence:", next_number)