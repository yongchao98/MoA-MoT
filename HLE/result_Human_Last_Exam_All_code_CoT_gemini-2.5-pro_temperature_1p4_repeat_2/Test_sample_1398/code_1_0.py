def solve():
    """
    Solves the grid transformation puzzle.
    The rule is determined to be that the output cell's value is the sum of all
    cells in the 3x3 neighborhood of the input cell (including the center), modulo 2.
    This is also known as the Parity rule in cellular automata.
    """
    grid = [
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ]

    rows = len(grid)
    cols = len(grid[0])
    output_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    def get_val(r, c):
        if 0 <= r < rows and 0 <= c < cols:
            return grid[r][c]
        return 0

    print("Equation: O(r,c) = ( Σ(i=r-1 to r+1) Σ(j=c-1 to c+1) I(i,j) ) mod 2\n")

    for r in range(rows):
        row_str = []
        for c in range(cols):
            # Sum of the 3x3 neighborhood
            total = 0
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    total += get_val(r + dr, c + dc)
            
            output_grid[r][c] = total % 2
            row_str.append(str(output_grid[r][c]))
        print(' '.join(row_str))

    # The problem asks for the output formatted as a single binary string.
    # We provide this as the final answer.
    binary_string = "".join(["".join(map(str, row)) for row in output_grid])
    print(f"\nFinal Answer String: {binary_string}")

solve()