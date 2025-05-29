def generate_grid():
    grid = [['❌' for _ in range(20)] for _ in range(20)]
    pattern = [
        [5, 11, 18],
        [3, 10, 17],
        [7, 14],
        [2, 9, 16],
        [6, 13],
        [1, 8, 15],
        [4, 11]
    ]
    
    for i in range(20):
        row_pattern = pattern[i % 7]
        for col in row_pattern:
            grid[i][col] = '✅'
    
    return grid

def print_grid(grid):
    for row in grid:
        print(''.join(row))

grid = generate_grid()
print_grid(grid)