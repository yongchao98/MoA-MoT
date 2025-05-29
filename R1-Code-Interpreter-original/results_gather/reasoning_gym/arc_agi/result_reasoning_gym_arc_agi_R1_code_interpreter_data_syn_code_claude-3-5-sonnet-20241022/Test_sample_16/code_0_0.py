def extract_pattern(input_grid):
    # Convert input string to 6x6 grid
    rows = input_grid.strip().split('\n')
    grid = [[int(x) for x in row.split()] for row in rows]
    
    # Create 3x3 output grid by taking elements from odd rows
    output = []
    for i in range(0, 6, 2):  # Take odd rows (0, 2, 4)
        row = [grid[i][0], grid[i][2], grid[i][4]]  # Take 1st, 3rd, 5th elements
        output.append(row)
    
    # Print the result in required format
    for row in output:
        print(' '.join(map(str, row)))

# Test with the new input
test_input = """9 6 3 6 3 6
6 6 6 6 6 6
9 6 6 6 6 6
6 6 6 6 6 6
1 6 0 6 0 6
6 6 6 6 6 6"""

extract_pattern(test_input)