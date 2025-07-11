def solve_sudoku(grid):
    """
    Solves a Sudoku puzzle using backtracking.
    """

    def find_empty(bo):
        """Finds an empty cell in the board."""
        for i in range(len(bo)):
            for j in range(len(bo[0])):
                if bo[i][j] == 0:
                    return (i, j)  # row, col
        return None

    def is_valid(bo, num, pos):
        """Checks if a number is valid in a given position."""
        # Check row
        for i in range(len(bo[0])):
            if bo[pos[0]][i] == num and pos[1] != i:
                return False

        # Check column
        for i in range(len(bo)):
            if bo[i][pos[1]] == num and pos[0] != i:
                return False

        # Check 3x3 box
        box_x = pos[1] // 3
        box_y = pos[0] // 3
        for i in range(box_y * 3, box_y * 3 + 3):
            for j in range(box_x * 3, box_x * 3 + 3):
                if bo[i][j] == num and (i, j) != pos:
                    return False
        return True

    def solve():
        """The recursive backtracking solver."""
        find = find_empty(grid)
        if not find:
            return True  # Puzzle is solved
        else:
            row, col = find

        for i in range(1, 10):
            if is_valid(grid, i, (row, col)):
                grid[row][col] = i

                if solve():
                    return True

                grid[row][col] = 0  # Backtrack

        return False

    solve()
    return grid

# Based on the visual representation in the problem description
initial_grid = [
    [5, 0, 0, 0, 0, 8, 4, 0, 9],
    [0, 0, 0, 5, 0, 0, 0, 3, 0],
    [0, 0, 6, 3, 7, 0, 0, 0, 1],
    [1, 0, 5, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 2, 0, 8, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 8],
    [7, 0, 0, 0, 0, 4, 1, 0, 5],
    [0, 0, 3, 0, 0, 2, 0, 0, 0],
    [4, 0, 9, 0, 0, 5, 0, 0, 3]
]

solved_grid = solve_sudoku(initial_grid)

# Print the top horizontal line of the solved sudoku
top_line = solved_grid[0]
print(' '.join(map(str, top_line)))