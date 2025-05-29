def create_output_grid(input_grid):
    # Convert input string to 3x3 grid
    input_rows = input_grid.strip().split('\n')
    grid_3x3 = [list(map(int, row.split())) for row in input_rows]
    
    # Create 6x6 output grid
    output = [[0 for _ in range(6)] for _ in range(6)]
    
    # Fill top-left quadrant (direct copy)
    for i in range(3):
        for j in range(3):
            output[i][j] = grid_3x3[i][j]
    
    # Fill top-right quadrant (horizontal mirror)
    for i in range(3):
        for j in range(3):
            output[i][5-j] = grid_3x3[i][j]
    
    # Fill bottom half (vertical mirror of top half)
    for i in range(3):
        for j in range(6):
            output[5-i][j] = output[i][j]
    
    # Convert to string format
    result = '\n'.join(' '.join(str(x) for x in row) for row in output)
    print(result)

# Test input
test_input = """9 9 5
7 7 4
7 7 7"""

create_output_grid(test_input)