def process_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]
    
    # First, process all 5s (down-right diagonals)
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 5:
                r, c = i, j
                while r < rows and c < cols:
                    if (r, c) != (i, j) and (input_grid[r][c] in [2, 5]):
                        break
                    output_grid[r][c] = 5
                    r += 1
                    c += 1
    
    # Then, process all 2s (up-right diagonals), overwriting 5s if needed
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 2:
                r, c = i, j
                while r >= 0 and c < cols:
                    if (r, c) != (i, j) and (input_grid[r][c] in [2, 5]):
                        break
                    output_grid[r][c] = 2
                    r -= 1
                    c += 1
    
    return output_grid

# Parse input grid
input_str = """4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 5 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 2 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 5 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 4
4 4 4 4 4 4 2 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 2 4 4 4 4 4 4 4 4 4 5 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4"""

# Convert input string to grid
input_grid = [[int(x) for x in line.split()] for line in input_str.split('\n')]

# Process grid and format output
output_grid = process_grid(input_grid)
output_str = '\n'.join(' '.join(str(x) for x in row) for row in output_grid)
print(output_str)