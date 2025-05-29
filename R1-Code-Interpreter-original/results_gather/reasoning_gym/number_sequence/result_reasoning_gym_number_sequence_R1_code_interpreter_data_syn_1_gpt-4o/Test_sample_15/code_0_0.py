sequence = [9, 4, 6, 5, 5, 5, 5]
differences = [sequence[i] - sequence[i+1] for i in range(len(sequence) - 1)]
next_number = sequence[-1]  # Assuming the sequence stabilizes at 5
print(next_number)