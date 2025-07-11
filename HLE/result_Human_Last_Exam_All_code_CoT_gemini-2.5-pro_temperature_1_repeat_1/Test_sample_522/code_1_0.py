def solve_puzzle(input_str):
    """
    Solves the puzzle by transforming the input grid according to the deduced rules.
    """
    # 1. Parse input into a grid
    rows_str = input_str.split(',')
    grid = [list(row) for row in rows_str]
    height = len(grid)
    if height == 0:
        print("MISSING")
        return
    width = len(grid[0])

    # 2. Find the coordinates (r2, c2) of '2'
    r2, c2 = -1, -1
    for r_idx, row in enumerate(grid):
        try:
            c_idx = row.index('2')
            r2, c2 = r_idx, c_idx
            break
        except ValueError:
            continue
    
    if r2 == -1: # '2' not found
        print("MISSING")
        return

    # 3. Determine direction based on parity of coordinates
    dr = -1 if r2 % 2 != 0 else 1
    dc = -1 if c2 % 2 != 0 else 1

    # 4. Trace path to find the first '0'
    cr, cc = r2, c2
    target_pos = None
    while True:
        cr += dr
        cc += dc

        # Check for out of bounds
        if not (0 <= cr < height and 0 <= cc < width):
            break

        # Check for target '0'
        if grid[cr][cc] == '0':
            target_pos = (cr, cc)
            break
    
    # 5. Determine the final output
    if target_pos is None:
        print("MISSING")
        return

    tr, tc = target_pos
    # Check for invalid boundary condition
    if tr == height - 1 or tc == width - 1:
        print("MISSING")
        return

    # 6. Perform the swap and print the result
    grid[r2][c2] = '0'
    grid[tr][tc] = '2'
    
    output_rows = ["".join(row) for row in grid]
    print(",".join(output_rows))


# The input from the puzzle description to find the value of "MISSING"
input_c = '000000,011120,111111'

# Solve for Input_C
solve_puzzle(input_c)