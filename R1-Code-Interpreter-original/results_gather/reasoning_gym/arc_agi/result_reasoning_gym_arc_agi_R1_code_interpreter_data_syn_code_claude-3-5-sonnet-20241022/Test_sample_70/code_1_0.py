def create_output_grid(input_grid):
    # Convert input string to 3x3 grid
    input_rows = input_grid.strip().split('\n')
    grid_3x3 = [list(map(int, row.split())) for row in input_rows]
    
    # Create 6x6 output grid
    output = [[0 for _ in range(6)] for _ in range(6)]
    
    # Fill first three rows and mirror horizontally
    for i in range(3):
        # Copy first three numbers
        for j in range(3):
            output[i][j] = grid_3x3[i][j]
        # Mirror horizontally
        for j in range(3):
            output[i][j+3] = grid_3x3[i][2-j]
    
    # Mirror vertically for bottom three rows
    for i in range(3):
        output[5-i] = output[i].copy()
    
    # Convert to string format
    result = '\n'.join(' '.join(str(x) for x in row) for row in output)
    print(result)

# Test input
test_input = """9 9 5
7 7 4
7 7 7"""

create_output_grid(test_input)