import numpy as np

def solve_grid(grid):
    rows, cols = grid.shape
    output_grid = np.zeros_like(grid)
    for r in range(rows):
        for c in range(cols):
            val = grid[r, c]
            same_color_neighbors = 0
            diff_color_neighbors = 0
            
            # Check orthogonal neighbors
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    neighbor_val = grid[nr, nc]
                    if neighbor_val == val:
                        same_color_neighbors += 1
                    else:
                        diff_color_neighbors += 1
            
            # Apply the transformation rule
            if diff_color_neighbors > same_color_neighbors:
                output_grid[r, c] = 1 - val
            else:
                output_grid[r, c] = val
    return output_grid

# Test Input Grid 6
input_grid_6 = np.array([
    [0, 1, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 0, 1, 0, 1],
    [0, 1, 0, 0, 1],
    [1, 0, 1, 1, 0]
])

output_grid_6 = solve_grid(input_grid_6)
result_string = "".join(map(str, output_grid_6.flatten()))
print(result_string)
