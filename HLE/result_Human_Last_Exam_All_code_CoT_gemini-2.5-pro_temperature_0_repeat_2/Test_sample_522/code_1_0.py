def solve_puzzle(input_str):
    """
    Solves the puzzle by transforming the grid based on the movement of the '2'.
    """
    # 1. Parse the input string into a 2D grid
    grid_rows = input_str.split(',')
    grid = [list(row) for row in grid_rows]
    
    # 2. Find the location of the '2'
    r2, c2 = -1, -1
    for r_idx, row in enumerate(grid):
        for c_idx, char in enumerate(row):
            if char == '2':
                r2, c2 = r_idx, c_idx
                break
        if r2 != -1:
            break
            
    # If no '2' is found, return the original string
    if r2 == -1:
        return input_str
        
    # 3. Determine the movement rule
    # The rule depends on the character to the left of the '2'
    # It's assumed '2' is not in the first column or first row based on examples
    v_left = grid[r2][c2 - 1]
    
    r_new = r2 - 1
    c_new = -1
    
    if v_left == '1':
        # Deflect left
        c_new = c2 - 1
    elif v_left == '0':
        # Deflect right
        c_new = c2 + 1
        
    # 4. Apply the transformation
    if c_new != -1:
        # Original position of '2' becomes '0'
        grid[r2][c2] = '0'
        # New position becomes '2'
        grid[r_new][c_new] = '2'
        
    # 5. Format the output
    output_rows = ["".join(row) for row in grid]
    result = ",".join(output_rows)
    
    return result

# Input for the puzzle to be solved
input_c = '000000,011120,111111'

# Calculate the value of MISSING
missing_value = solve_puzzle(input_c)

# Print the final result
print(missing_value)