def process_grid(input_grid_str):
    # Parse input grid
    input_rows = input_grid_str.strip().split('\n')
    grid = [list(map(int, row.split())) for row in input_rows]
    
    # Create output by duplicating the input grid
    output_grid = grid + grid  # Concatenate the grid with itself
    
    # Format output
    output_str = ''
    for row in output_grid:
        output_str += ' '.join(map(str, row)) + '\n'
    
    print(output_str.strip())

# Test input
test_input = """7 6 6 6 6
6 7 1 6 6
6 1 9 9 6
6 1 6 8 8"""

process_grid(test_input)