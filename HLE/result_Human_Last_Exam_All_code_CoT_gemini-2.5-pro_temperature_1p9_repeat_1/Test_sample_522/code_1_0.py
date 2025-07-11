def solve_puzzle():
    """
    Solves the puzzle by transforming the input grid based on a discovered rule.
    The rule is that a '2' moves diagonally up-right if the target cell is '0',
    otherwise it tries to move diagonally up-left if that target cell is '0'.
    """
    input_c = '000000,011120,111111'

    # Step 1: Parse the input string into a 2D grid
    rows_str = input_c.split(',')
    grid = [list(row) for row in rows_str]

    if not grid or not grid[0]:
        print(input_c)
        return

    num_rows = len(grid)
    num_cols = len(grid[0])

    # Step 2: Find the location of '2'. Assuming one '2' per puzzle.
    r, c = -1, -1
    found_2 = False
    for i in range(num_rows):
        for j in range(num_cols):
            if grid[i][j] == '2':
                r, c = i, j
                found_2 = True
                break
        if found_2:
            break

    # If no '2' is found, the grid doesn't change
    if not found_2:
        print(input_c)
        return
        
    # Step 3: Apply the move logic
    moved = False
    
    # Priority 1: Check up-right move
    # Check bounds: r > 0 and c+1 is within the grid width
    if r > 0 and c + 1 < num_cols:
        if grid[r - 1][c + 1] == '0':
            grid[r - 1][c + 1] = '2'
            grid[r][c] = '0'
            moved = True
            
    # Priority 2: If not moved, check up-left move
    if not moved:
        # Check bounds: r > 0 and c > 0
        if r > 0 and c > 0:
            if grid[r - 1][c - 1] == '0':
                grid[r - 1][c - 1] = '2'
                grid[r][c] = '0'
                moved = True
                
    # Step 4: Format the resulting grid back into a string and print
    output_rows = ["".join(row) for row in grid]
    result = ",".join(output_rows)
    print(result)

solve_puzzle()