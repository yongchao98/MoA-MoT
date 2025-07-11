def solve_puzzle():
    """
    Solves the puzzle by transforming the input grid based on the movement of the digit '2'.
    """
    input_str = '000000,011120,111111'

    # Parse the input string into a 2D grid of characters
    rows = input_str.split(',')
    grid = [list(row) for row in rows]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0

    # Find the location of '2'
    pos_2 = None
    for r_idx, row in enumerate(grid):
        try:
            c_idx = row.index('2')
            pos_2 = (r_idx, c_idx)
            break
        except ValueError:
            continue

    # If '2' is not found, there is nothing to do.
    if pos_2 is None:
        print(input_str)
        return
        
    r, c = pos_2
    
    # Rule 1: Attempt to move diagonally up-right
    # Check if target (r-1, c+1) is within bounds and is '0'
    if r > 0 and c + 1 < width and grid[r - 1][c + 1] == '0':
        grid[r][c] = '0'
        grid[r - 1][c + 1] = '2'
    # Rule 2: Attempt to move diagonally up-left
    # Check if target (r-1, c-1) is within bounds and is '0'
    elif r > 0 and c > 0 and grid[r - 1][c - 1] == '0':
        grid[r][c] = '0'
        grid[r - 1][c - 1] = '2'

    # Convert the modified grid back to a string and print it
    output_rows = [''.join(row) for row in grid]
    result = ','.join(output_rows)
    print(result)

solve_puzzle()