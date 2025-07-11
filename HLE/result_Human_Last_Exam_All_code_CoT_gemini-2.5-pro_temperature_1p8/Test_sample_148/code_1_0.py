def solve_sudoku(grid):
    """
    Solves a Sudoku puzzle using backtracking.
    """
    find = find_empty(grid)
    if not find:
        return True  # Puzzle is solved
    else:
        row, col = find

    for num in range(1, 10):
        if is_valid(grid, row, col, num):
            grid[row][col] = num

            if solve_sudoku(grid):
                return True

            # Backtrack if the solution is not found
            grid[row][col] = 0

    return False

def find_empty(grid):
    """
    Finds the first empty cell (represented by 0) in the grid.
    """
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return (i, j)  # row, col
    return None

def is_valid(grid, row, col, num):
    """
    Checks if placing a number in a given cell is valid.
    """
    # Check if the number is already in the row
    for x in range(9):
        if grid[row][x] == num:
            return False

    # Check if the number is already in the column
    for x in range(9):
        if grid[x][col] == num:
            return False

    # Check if the number is in the 3x3 sub-grid
    start_row = row - row % 3
    start_col = col - col % 3
    for i in range(3):
        for j in range(3):
            if grid[i + start_row][j + start_col] == num:
                return False
                
    return True

if __name__ == '__main__':
    # The Sudoku grid parsed from the input
    puzzle = [
        [5, 0, 0, 0, 8, 0, 0, 4, 9],
        [0, 0, 0, 5, 0, 0, 0, 3, 0],
        [0, 6, 7, 3, 0, 0, 0, 0, 1],
        [1, 0, 5, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 8, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 8],
        [7, 0, 0, 0, 0, 4, 1, 5, 0],
        [0, 3, 0, 0, 0, 2, 0, 0, 0],
        [4, 0, 9, 0, 5, 0, 0, 0, 3]
    ]

    if solve_sudoku(puzzle):
        # Print the top horizontal line of the solved sudoku
        # The instruction asks to output each number in the final equation.
        # This is interpreted as printing the numbers of the solved first row.
        print(' '.join(map(str, puzzle[0])))
    else:
        print("No solution exists for the given Sudoku puzzle.")