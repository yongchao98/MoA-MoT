# Define the grid
grid = [
    [1, 2],
    [3, 2]
]

# Pre-selected position
selected_positions = [(0, 1)]

# Check possible selections
for i in range(2):
    for j in range(2):
        if (i, j) not in selected_positions:
            selected_positions.append((i, j))
            if len(selected_positions) == 3:
                # Calculate row and column sums
                row_sums = [0, 0]
                col_sums = [0, 0]
                for pos in selected_positions:
                    row_sums[pos[0]] += grid[pos[0]][pos[1]]
                    col_sums[pos[1]] += grid[pos[0]][pos[1]]
                
                # Check if all sums are <= 4
                if all(sum_val <= 4 for sum_val in row_sums + col_sums):
                    print(selected_positions)
                    break
            selected_positions.pop()