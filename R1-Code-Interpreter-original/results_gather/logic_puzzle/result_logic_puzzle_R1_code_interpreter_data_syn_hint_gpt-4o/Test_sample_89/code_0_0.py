grid = [
    [8, 12, 9, 3, -11, 9],
    [12, -12, -4, 11, 10, 7],
    [-8, 6, 9, 10, 1, 3],
    [5, 11, 6, 7, -9, 11],
    [6, 1, 11, 12, 4, -3],
    [2, 9, 4, 3, 10, 9]
]

selected_positions = [(0, 5)]  # Starting with the given position

# We need to select 11 more numbers
# Strategy: Select positive numbers first, then ensure even negatives

# Track the number of negatives in each row and column
row_negatives = [0] * 6
col_negatives = [0] * 6

# Mark the initial selection
row_negatives[0] += 0 if grid[0][5] >= 0 else 1
col_negatives[5] += 0 if grid[0][5] >= 0 else 1

# Function to check if we can select a number
def can_select(r, c):
    if (r, c) in selected_positions:
        return False
    if grid[r][c] < 0:
        # Check if adding this negative number keeps the negatives even
        if (row_negatives[r] + 1) % 2 != 0 or (col_negatives[c] + 1) % 2 != 0:
            return False
    return True

# Select numbers
for r in range(6):
    for c in range(6):
        if len(selected_positions) >= 12:
            break
        if can_select(r, c):
            selected_positions.append((r, c))
            if grid[r][c] < 0:
                row_negatives[r] += 1
                col_negatives[c] += 1

# Output the selected positions
print(', '.join(f"{r} {c}" for r, c in selected_positions))