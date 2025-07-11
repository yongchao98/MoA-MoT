def solve_puzzle(input_str):
    """
    Solves the puzzle by transforming the grid based on the position of '2'.

    The rule is:
    1. A '2' at position (r, c) becomes '0'.
    2. A '0' at a target diagonal position becomes '2'.
    3. The direction of the diagonal jump depends on the parity of (r, c):
       - (odd, odd):   up-left (-1, -1)
       - (odd, even):  up-right (-1, +1)
       - (even, odd):  down-left (+1, -1)
       - (even, even): down-right (+1, +1)
    4. The transformation only happens if the target cell is within bounds and contains a '0'.
    """
    # 1. Parse the input string into a 2D list of integers
    rows_str = input_str.split(',')
    grid = [[int(char) for char in row] for row in rows_str]
    
    # Create a copy to modify, preserving the original for lookups
    output_grid = [row[:] for row in grid]
    
    num_rows = len(grid)
    if num_rows == 0:
        print(input_str)
        return
    num_cols = len(grid[0])

    # 2. Find all locations of the number 2
    locs_of_2 = []
    for r in range(num_rows):
        for c in range(num_cols):
            if grid[r][c] == 2:
                locs_of_2.append((r, c))

    # 3. Apply the transformation for each '2' found
    for r, c in locs_of_2:
        dr, dc = 0, 0
        # Determine jump vector based on parity of coordinates
        if r % 2 != 0:  # Odd row
            if c % 2 != 0:  # Odd col
                dr, dc = -1, -1
            else:  # Even col
                dr, dc = -1, 1
        else:  # Even row
            if c % 2 != 0:  # Odd col
                dr, dc = 1, -1
            else:  # Even col
                dr, dc = 1, 1
        
        # Calculate target coordinates
        tr, tc = r + dr, c + dc

        # Check if target is valid and contains a 0
        if 0 <= tr < num_rows and 0 <= tc < num_cols and grid[tr][tc] == 0:
            # Apply the swap
            output_grid[r][c] = 0
            output_grid[tr][tc] = 2

    # 4. Convert the output grid back to a string and print
    output_rows = [''.join(map(str, row)) for row in output_grid]
    result = ','.join(output_rows)
    print(result)

# The input for which we need to find the 'MISSING' value
input_c = '000000,011120,111111'
solve_puzzle(input_c)