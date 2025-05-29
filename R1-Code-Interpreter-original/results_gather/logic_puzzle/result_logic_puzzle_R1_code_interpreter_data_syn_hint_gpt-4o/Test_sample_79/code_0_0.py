grid = [
    [-1, 2, -3, 8, -9],
    [2, -8, -5, 1, -3],
    [-10, 6, 1, 1, 10],
    [-2, -2, 6, 9, -6],
    [-3, -10, 2, 6, -9]
]

selected_positions = [(0, 4)]
selected_numbers = [-9]

# We need to select 9 more numbers
for i in range(5):
    for j in range(5):
        if len(selected_positions) == 10:
            break
        if (i, j) not in selected_positions:
            # Check if adding this number keeps the product of each row and column positive
            row_neg_count = sum(1 for x in range(5) if (i, x) in selected_positions and grid[i][x] < 0)
            col_neg_count = sum(1 for x in range(5) if (x, j) in selected_positions and grid[x][j] < 0)
            
            if grid[i][j] < 0:
                row_neg_count += 1
                col_neg_count += 1
            
            if row_neg_count % 2 == 0 and col_neg_count % 2 == 0:
                selected_positions.append((i, j))
                selected_numbers.append(grid[i][j])

# Output the selected positions
print(selected_positions)