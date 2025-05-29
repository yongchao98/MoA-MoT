# List of column positions of check marks
positions = [17, 18, 12, 9, 4, 11, 6, 16, 14, 5, 15, 2, 3, 13, 10, 7, 1, 8, 19]

# Find the missing number in the sequence from 1 to 20
missing_number = set(range(1, 21)) - set(positions)
print(missing_number)