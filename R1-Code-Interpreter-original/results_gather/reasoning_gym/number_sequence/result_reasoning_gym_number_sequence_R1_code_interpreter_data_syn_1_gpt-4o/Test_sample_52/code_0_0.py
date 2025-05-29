# Given sequence
sequence = [4, 15, 30, 56, 97]

# Calculate first-level differences
first_differences = [sequence[i] - sequence[i - 1] for i in range(1, len(sequence))]

# Calculate second-level differences
second_differences = [first_differences[i] - first_differences[i - 1] for i in range(1, len(first_differences))]

# Assume the pattern in second-level differences continues
next_second_difference = second_differences[-1] + (second_differences[-1] - second_differences[-2])

# Calculate the next first-level difference
next_first_difference = first_differences[-1] + next_second_difference

# Calculate the next term in the sequence
next_term = sequence[-1] + next_first_difference

print(next_term)