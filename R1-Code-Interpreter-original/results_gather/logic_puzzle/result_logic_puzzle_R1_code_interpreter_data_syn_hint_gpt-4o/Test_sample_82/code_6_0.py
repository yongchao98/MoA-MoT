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

# Initialize row and column negative counts
row_neg_count = [0] * 6
col_neg_count = [0] * 6

# Mark the initial selection
row_neg_count[0] += 0
col_neg_count[5] += 0

# Function to check if adding a number keeps the product positive
def is_valid_selection(row_neg_count, col_neg_count, row, col, num):
    if num < 0:
        row_neg_count[row] += 1
        col_neg_count[col] += 1
    valid = (row_neg_count[row] % 2 == 0) and (col_neg_count[col] % 2 == 0)
    if num < 0:
        row_neg_count[row] -= 1
        col_neg_count[col] -= 1
    return valid

# Backtracking function to find a valid selection
def find_selection(count):
    if count == 12:
        return True
    for i in range(6):
        for j in range(6):
            if (i, j) not in selected_positions:
                num = grid[i][j]
                if num > 0 or is_valid_selection(row_neg_count, col_neg_count, i, j, num):
                    selected_positions.append((i, j))
                    selected_numbers.append(num)
                    if num < 0:
                        row_neg_count[i] += 1
                        col_neg_count[j] += 1
                    if find_selection(count + 1):
                        return True
                    selected_positions.pop()
                    selected_numbers.pop()
                    if num < 0:
                        row_neg_count[i] -= 1
                        col_neg_count[j] -= 1
    return False

# Start the selection process
find_selection(1)

# Output the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))