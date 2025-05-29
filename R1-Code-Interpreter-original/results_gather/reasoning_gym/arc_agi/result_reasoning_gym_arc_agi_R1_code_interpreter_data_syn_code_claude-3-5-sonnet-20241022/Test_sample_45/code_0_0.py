def find_bounds(grid):
    rows, cols = len(grid), len(grid[0])
    min_row = min_col = float('inf')
    max_row = max_col = -1
    
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] != 0:
                min_row = min(min_row, i)
                max_row = max(max_row, i)
                min_col = min(min_col, j)
                max_col = max(max_col, j)
    
    return min_row, max_row, min_col, max_col

def process_grid(input_grid):
    # Find bounds of non-zero numbers
    min_row, max_row, min_col, max_col = find_bounds(input_grid)
    
    # Extract the relevant subgrid
    output = []
    for i in range(min_row, max_row + 1):
        row = []
        for j in range(min_col, max_col + 1):
            if input_grid[i][j] != 0:
                row.append(input_grid[i][j])
            else:
                # Fill gaps with average of non-zero neighbors
                neighbors = []
                for di in [-1, 0, 1]:
                    for dj in [-1, 0, 1]:
                        ni, nj = i + di, j + dj
                        if (0 <= ni < len(input_grid) and 
                            0 <= nj < len(input_grid[0]) and 
                            input_grid[ni][nj] != 0):
                            neighbors.append(input_grid[ni][nj])
                row.append(round(sum(neighbors)/len(neighbors)) if neighbors else 0)
        output.append(row)
    
    # Print the output grid
    for row in output:
        print(*row)

# Test input
test_input = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 5, 2, 0, 0],
    [0, 0, 2, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
]

process_grid(test_input)