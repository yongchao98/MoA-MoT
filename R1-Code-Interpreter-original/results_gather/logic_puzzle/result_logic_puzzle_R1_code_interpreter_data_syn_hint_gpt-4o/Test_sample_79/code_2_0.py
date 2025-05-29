grid = [
    [-1, 2, -3, 8, -9],
    [2, -8, -5, 1, -3],
    [-10, 6, 1, 1, 10],
    [-2, -2, 6, 9, -6],
    [-3, -10, 2, 6, -9]
]

selected_positions = [(0, 4)]
selected_numbers = [-9]

# Track negative counts in rows and columns
row_neg_count = [0] * 5
col_neg_count = [0] * 5
row_neg_count[0] += 1
col_neg_count[4] += 1

# Function to check if adding a number keeps the product condition valid
def can_add(i, j):
    temp_row_neg = row_neg_count[i]
    temp_col_neg = col_neg_count[j]
    if grid[i][j] < 0:
        temp_row_neg += 1
        temp_col_neg += 1
    return temp_row_neg % 2 == 0 and temp_col_neg % 2 == 0

# Greedily select numbers
for i in range(5):
    for j in range(5):
        if len(selected_positions) == 10:
            break
        if (i, j) not in selected_positions and can_add(i, j):
            selected_positions.append((i, j))
            selected_numbers.append(grid[i][j])
            if grid[i][j] < 0:
                row_neg_count[i] += 1
                col_neg_count[j] += 1

print(selected_positions)