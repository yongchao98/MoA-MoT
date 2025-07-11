def solve_sudoku():
    """
    Solves the Sudoku puzzle provided in the problem description.
    The puzzle is first parsed from its visual format into a 9x9 grid.
    Then, a backtracking algorithm finds the solution.
    Finally, it prints the top horizontal line of the solved puzzle.
    """
    # The grid is parsed from the LaTeX-style matrix in the prompt.
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
        """Finds an empty cell in the board"""
        for i in range(len(bo)):
            for j in range(len(bo[0])):
                if bo[i][j] == 0:
                    return (i, j)  # row, col
        return None

    def is_valid(bo, num, pos):
        """Checks if a number is valid in a given position"""
        # Check row
        for i in range(len(bo[0])):
            if bo[pos[0]][i] == num and pos[1] != i:
                return False

        # Check column
        for i in range(len(bo)):
            if bo[i][pos[1]] == num and pos[0] != i:
                return False

        # Check box
        box_x = pos[1] // 3
        box_y = pos[0] // 3

        for i in range(box_y*3, box_y*3 + 3):
            for j in range(box_x * 3, box_x*3 + 3):
                if bo[i][j] == num and (i,j) != pos:
                    return False
        return True

    def solve():
        """
        Solves the sudoku board using backtracking
        """
        find = find_empty(grid)
        if not find:
            return True
        else:
            row, col = find

        for i in range(1,10):
            if is_valid(grid, i, (row, col)):
                grid[row][col] = i

                if solve():
                    return True

                grid[row][col] = 0
        return False

    solve()
    
    # Print the top horizontal line of the solved sudoku
    # The problem asks to output each number in the final equation.
    # This is interpreted as printing the space-separated numbers of the first row.
    print(' '.join(map(str, grid[0])))

solve_sudoku()