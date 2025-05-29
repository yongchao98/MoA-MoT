# Generate the sequence to verify the pattern
sequence = [1, 0, 0, 0]
# Extend the sequence to check the pattern
for i in range(4, 8):
    sequence.append(sequence[i % 4])

# Output the fifth element in the sequence
print(sequence[4])