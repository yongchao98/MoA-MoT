def analyze_grid(input_grid):
    # Convert input string to 2D grid
    rows = input_grid.strip().split('\n')
    grid = [[int(x) for x in row.split()] for row in rows]
    rows, cols = len(grid), len(grid[0])
    
    # Divide grid into 9 regions and find dominant non-3 numbers
    height, width = rows//3, cols//3
    output = [[3]*3 for _ in range(3)]
    
    def get_region_special(region_rows, region_cols):
        special_counts = {}
        for r in range(region_rows[0], region_rows[1]):
            for c in range(region_cols[0], region_cols[1]):
                if grid[r][c] != 3:
                    special_counts[grid[r][c]] = special_counts.get(grid[r][c], 0) + 1
        if special_counts:
            # Return the most frequent special number
            return max(special_counts.items(), key=lambda x: x[1])[0]
        return 3
    
    # Process each region
    for i in range(3):
        for j in range(3):
            row_range = (i*height, (i+1)*height)
            col_range = (j*width, (j+1)*width)
            special = get_region_special(row_range, col_range)
            output[i][j] = special
    
    # Apply symmetry rules based on the examples
    # If we have 4s in middle, they should spread to corners
    if output[1][1] == 4:
        output[0][0] = output[0][2] = output[2][0] = output[2][2] = 2
    
    # Print the output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = """3 3 3 3 3 2 3 3 3 3 2 2 3 3 2 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 2 3 3 3 3 3 3 3 3
2 3 3 3 2 3 3 3 3 2 3 3 2 3 3 3 3 3 3 3 3 2 3 3 3 3
3 3 3 3 3 3 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 2 3 3 3 3 2 3 2 3 3 3 3 2 2 3 2 3
2 3 3 3 3 3 2 3 2 3 3 3 4 4 4 4 3 3 3 3 2 3 3 2 3 3
3 3 3 3 3 3 2 3 3 3 3 3 4 4 4 4 3 2 3 3 3 3 3 3 3 2
2 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 3 3 3 2 3 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 2 3 3 3 3 2 3 3 2 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 3 3 2 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3
3 2 3 3 3 3 3 3 3 2 2 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 3 3 3 3 3 3
3 3 2 3 2 3 3 3 3 3 3 2 4 4 4 4 3 3 3 3 4 4 4 4 3 3
3 3 3 2 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 4 4 4 4 3 3
3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 3 3 3 3 4 4 4 4 3 3
3 3 3 3 3 3 2 3 3 3 3 3 4 4 4 4 3 3 3 3 4 4 4 4 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 2 3 3 3 3 2 3 3 3 3
3 3 3 3 2 2 2 2 3 3 3 3 2 3 3 3 3 2 3 3 3 3 3 3 3 3
3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 3 3 3 3 2 3 3 3
3 3 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3"""

analyze_grid(test_input)