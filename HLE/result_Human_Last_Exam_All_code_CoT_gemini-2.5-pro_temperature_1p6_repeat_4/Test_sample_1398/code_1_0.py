def solve():
    """
    Applies the transformation rule to the input grid.
    The rule is as follows:
    1. For each cell, find its 4 diagonal neighbors (padding with 0 outside the grid).
    2. Sum the values of these neighbors to get S_diag.
    3. If the input cell value is 0, the output cell is 1.
    4. If the input cell value is 1, the output is 1 if S_diag is even, and 0 if S_diag is odd.
    """
    grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    R, C = len(grid), len(grid[0])
    out_grid = [[0 for _ in range(C)] for _ in range(R)]

    for r in range(R):
        for c in range(C):
            s_diag = 0
            # NW
            if r > 0 and c > 0: s_diag += grid[r-1][c-1]
            # NE
            if r > 0 and c < C - 1: s_diag += grid[r-1][c+1]
            # SW
            if r < R - 1 and c > 0: s_diag += grid[r+1][c-1]
            # SE
            if r < R - 1 and c < C - 1: s_diag += grid[r+1][c+1]

            if grid[r][c] == 0:
                out_grid[r][c] = 1
            else:  # grid[r][c] == 1
                if s_diag % 2 == 0:
                    out_grid[r][c] = 1
                else:
                    out_grid[r][c] = 0

    # Print the resulting grid as a flat binary string
    output_string = "".join(str(cell) for row in out_grid for cell in row)
    print(output_string)

solve()