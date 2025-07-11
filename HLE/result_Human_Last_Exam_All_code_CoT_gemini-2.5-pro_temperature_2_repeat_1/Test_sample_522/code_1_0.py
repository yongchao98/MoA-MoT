def solve_puzzle(input_str):
    """
    Solves the puzzle based on the derived rules.
    The primary rule is that if a '2' is found at column index W-2 (where W is grid width),
    the output is 'MISSING'. Otherwise, a transformation would occur.
    This function specifically calculates the result for the provided inputs.
    """
    rows_str = input_str.split(',')
    grid = [list(row) for row in rows_str]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0

    pos_2 = None
    for r in range(height):
        for c in range(width):
            if grid[r][c] == '2':
                pos_2 = (r, c)
                break
        if pos_2:
            break

    # If no '2' is found, no transformation happens.
    if not pos_2:
        print(input_str)
        return

    r, c = pos_2

    # Rule: If '2' is at column index (width - 2), output is "MISSING".
    if c == width - 2:
        print("MISSING")
        return

    # --- Full transformation logic for cases A and B (not strictly needed for C) ---
    # Priority: TR > TL
    moved = False
    
    # Try TR move
    if r > 0 and c < width - 1:
        tr_r, tr_c = r - 1, c + 1
        t_r, t_c = r - 1, c
        r_r, r_c = r, c + 1
        if grid[tr_r][tr_c] == '0':
            t_val = grid[t_r][t_c] if t_r >= 0 else ''
            r_val = grid[r_r][r_c] if r_c < width else ''
            if t_val == '1' or r_val == '1':
                grid[r][c] = '0'
                grid[tr_r][tr_c] = '2'
                moved = True

    # Try TL move if TR didn't happen
    if not moved and r > 0 and c > 0:
        tl_r, tl_c = r - 1, c - 1
        t_r, t_c = r - 1, c
        l_r, l_c = r, c - 1
        if grid[tl_r][tl_c] == '0':
            t_val = grid[t_r][t_c] if t_r >= 0 else ''
            l_val = grid[l_r][l_c] if l_c >= 0 else ''
            if t_val == '1' or l_val == '1':
                grid[r][c] = '0'
                grid[tl_r][tl_c] = '2'
                moved = True
                
    if moved:
        result = ",".join(["".join(row) for row in grid])
        print(result)
    else:
        # If no valid move was found after all checks.
        print("MISSING")


# Calculate the value for Input_C
input_c = '000000,011120,111111'
solve_puzzle(input_c)