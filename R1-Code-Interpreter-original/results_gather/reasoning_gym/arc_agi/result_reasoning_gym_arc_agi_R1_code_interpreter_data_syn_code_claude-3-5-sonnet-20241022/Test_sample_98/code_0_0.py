def process_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Create a copy
    
    # Find all 5s and 2s
    fives = []
    twos = []
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 5:
                fives.append((i, j))
            elif input_grid[i][j] == 2:
                twos.append((i, j))
    
    # Process 5s (down-right diagonals)
    for five_pos in fives:
        i, j = five_pos
        r, c = i, j
        while r < rows and c < cols:
            # Stop if we hit a 2 or another 5
            if (r, c) != (i, j) and (input_grid[r][c] == 2 or input_grid[r][c] == 5):
                break
            output_grid[r][c] = 5
            r += 1
            c += 1
    
    # Process 2s (up-right diagonals)
    for two_pos in twos:
        i, j = two_pos
        r, c = i, j
        while r >= 0 and c < cols:
            # Stop if we hit a 5 or another 2
            if (r, c) != (i, j) and (input_grid[r][c] == 5 or input_grid[r][c] == 2):
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