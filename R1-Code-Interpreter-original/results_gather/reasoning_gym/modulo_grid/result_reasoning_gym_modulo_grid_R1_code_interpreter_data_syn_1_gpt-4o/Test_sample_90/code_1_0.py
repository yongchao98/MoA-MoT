# Define the column positions of the checkmarks
column_positions = [7, 12, 8, 6, 15, 4, 1, 18, 3, 14, 16, 13, 2, 19, 11, 9, 5, 10, 7, 13]

# Check for a pattern in the column positions
def find_column_pattern(columns):
    differences = []
    for i in range(1, len(columns)):
        differences.append(columns[i] - columns[i-1])
    return differences

column_differences = find_column_pattern(column_positions)
print("Column differences:", column_differences)