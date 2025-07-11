def solve_puzzle(input_str):
    """
    Solves the puzzle by finding a '2' and swapping it with a '0' neighbor
    based on a conditional priority rule.
    """
    # Parse the input string into a 2D list of integers
    grid = [list(map(int, row)) for row in input_str.split(',')]
    rows = len(grid)
    cols = len(grid[0])

    # Find the coordinates (r, c) of the number 2
    r2, c2 = -1, -1
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == 2:
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    # If no '2' is found, return the original grid (or handle as an error)
    if r2 == -1:
        return input_str

    # Define the two priority lists for directions (dr, dc)
    # P1: Up-Left, Up, Up-Right, Left, Right, Down-Left, Down, Down-Right
    priority1 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    # P2: Up, Up-Left, Up-Right, ...
    priority2 = [(-1, 0), (-1, -1), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

    # Check if r == c to select the appropriate priority list
    if r2 == c2:
        directions = priority1
    else:
        directions = priority2

    # Find the target '0' neighbor and perform the swap
    for dr, dc in directions:
        nr, nc = r2 + dr, c2 + dc

        # Check if the neighbor is within grid boundaries
        if 0 <= nr < rows and 0 <= nc < cols:
            # Check if the neighbor is a '0'
            if grid[nr][nc] == 0:
                # Perform the swap
                grid[r2][c2], grid[nr][nc] = grid[nr][nc], grid[r2][c2]
                break
    
    # Convert the modified grid back to the output string format
    output_rows = [ "".join(map(str, row)) for row in grid]
    print(','.join(output_rows))

# Solve for the MISSING value from Input_C
solve_puzzle('000000,011120,111111')