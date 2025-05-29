# Define the sequence
sequence = [5, 4, 4, 4, 4, 4]

# Check the pattern and predict the next number
if sequence[1:] == [4] * (len(sequence) - 1):
    next_number = 4
else:
    next_number = None

print(next_number)