# Given sequence
sequence = [3, 2, 4, 5]

# Calculate differences between consecutive numbers
differences = [sequence[i+1] - sequence[i] for i in range(len(sequence) - 1)]

# Print the differences
print("Differences:", differences)

# Predict the next number in the sequence
# Assuming a pattern in differences, let's try to find the next difference
# and compute the next number in the sequence

# The next difference should be the next in the observed pattern
next_difference = differences[len(differences) % len(differences)]

# Calculate the next number
next_number = sequence[-1] + next_difference

# Output the next number
print("Next number in the sequence:", next_number)