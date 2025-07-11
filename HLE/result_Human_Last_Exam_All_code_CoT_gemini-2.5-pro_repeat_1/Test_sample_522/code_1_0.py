def solve_puzzle():
    """
    Solves the puzzle by transforming the input grid based on a derived rule.
    The rule is that the '2' moves to a neighboring '0' based on the
    parity of its current coordinates.
    """
    input_str = '000000,011120,111111'

    # 1. Parse the input string into a 2D list of characters
    rows = input_str.split(',')
    grid = [list(row) for row in rows]

    # 2. Find the coordinates (r, c) of the '2'
    r2, c2 = -1, -1
    for r, row_list in enumerate(grid):
        if '2' in row_list:
            r2 = r
            c2 = row_list.index('2')
            break
    
    # This check is for robustness, assuming '2' is always present per puzzle rules.
    if r2 == -1:
        print("Error: '2' not found in the input grid.")
        return

    # 3. Determine the move vector (dr, dc) based on the parity of the coordinates.
    dr, dc = 0, 0
    is_r_odd = r2 % 2 != 0
    is_c_odd = c2 % 2 != 0

    if is_r_odd and is_c_odd:         # (odd, odd)
        dr, dc = -1, -1
    elif is_r_odd and not is_c_odd:   # (odd, even)
        dr, dc = -1, 1
    elif not is_r_odd and is_c_odd:   # (even, odd) - rule hypothesized
        dr, dc = 1, 1
    elif not is_r_odd and not is_c_odd: # (even, even) - rule hypothesized
        dr, dc = 1, -1

    # 4. Calculate the new coordinates for the '2'
    new_r, new_c = r2 + dr, c2 + dc

    # 5. Update the grid: move the '2' and replace its old spot with '0'
    # Assuming the target cell is always a valid '0' as per the puzzle's nature.
    if 0 <= new_r < len(grid) and 0 <= new_c < len(grid[0])):
        grid[r2][c2] = '0'
        grid[new_r][new_c] = '2'

    # 6. Format the grid back into a comma-separated string and print it
    output_rows = ["".join(row) for row in grid]
    result = ",".join(output_rows)
    print(result)

solve_puzzle()