grid = [
    [9, 5, 1, 12, -2, 11],
    [-5, 6, 6, 11, 10, 8],
    [8, -1, -6, 1, -8, 5],
    [1, 12, 11, 10, -9, 12],
    [3, 3, -12, -4, 5, 12],
    [-1, -10, -3, -2, 4, 12]
]

# Start with the given number
selected_positions = [(0, 5)]
selected_numbers = [11]

# Function to check if adding a number keeps the product positive
def is_valid_selection(row_counts, col_counts, row, col, num):
    if num < 0:
        row_counts[row] += 1
        col_counts[col] += 1
    valid = (row_counts[row] % 2 == 0) and (col_counts[col] % 2 == 0)
    if num < 0:
        row_counts[row] -= 1
        col_counts[col] -= 1
    return valid

# Initialize row and column negative counts
row_counts = [0] * 6
col_counts = [0] * 6

# Mark the initial selection
row_counts[0] += 0
col_counts[5] += 0

# Select numbers
for r in range(6):
    for c in range(6):
        if len(selected_positions) == 12:
            break
        if (r, c) not in selected_positions:
            num = grid[r][c]
            if num > 0 or is_valid_selection(row_counts, col_counts, r, c, num):
                selected_positions.append((r, c))
                selected_numbers.append(num)
                if num < 0:
                    row_counts[r] += 1
                    col_counts[c] += 1

# Output the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))