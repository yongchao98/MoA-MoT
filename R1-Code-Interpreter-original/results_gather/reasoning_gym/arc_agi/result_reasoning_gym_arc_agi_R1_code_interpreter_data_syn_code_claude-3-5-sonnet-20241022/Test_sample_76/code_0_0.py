def analyze_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Convert input string to 2D grid
    grid = []
    numbers = [int(x) for x in input_grid.split()]
    for i in range(0, len(numbers), cols):
        grid.append(numbers[i:i+cols])
    
    # Divide the grid into 9 regions (3x3 output grid)
    height = rows // 3
    width = cols // 3
    
    # Calculate the most frequent non-3 number for each region
    output = []
    for i in range(3):
        row = []
        for j in range(3):
            region = []
            for r in range(i*height, (i+1)*height):
                for c in range(j*width, (j+1)*width):
                    if r < rows and c < cols:
                        region.append(grid[r][c])
            
            # Find most common non-3 number
            non_3 = [x for x in region if x != 3]
            if non_3:
                most_common = max(set(non_3), key=non_3.count)
                row.append(most_common)
            else:
                row.append(3)
        output.append(row)
    
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