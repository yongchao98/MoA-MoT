def solve_sudoku(grid):
    """
    Solves a Sudoku puzzle using a backtracking algorithm.
    The input grid is modified in-place to contain the solution.
    """
    find = find_empty(grid)
    if not find:
        return True  # Puzzle is solved
    else:
        row, col = find

    for num in range(1, 10):
        if is_valid(grid, num, (row, col)):
            grid[row][col] = num

            if solve_sudoku(grid):
                return True

            grid[row][col] = 0  # Backtrack

    return False

def find_empty(grid):
    """
    Finds the first empty cell (represented by 0) in the grid.
    Returns a (row, col) tuple or None if no empty cell is found.
    """
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                return (i, j)
    return None

def is_valid(grid, num, pos):
    """
    Checks if placing a number in a given position is valid according
    to Sudoku rules.
    """
    row, col = pos

    # Check row
    for i in range(9):
        if grid[row][i] == num and col != i:
            return False

    # Check column
    for i in range(9):
        if grid[i][col] == num and row != i:
            return False

    # Check 3x3 box
    box_x = col // 3
    box_y = row // 3

    for i in range(box_y * 3, box_y * 3 + 3):
        for j in range(box_x * 3, box_x * 3 + 3):
            if grid[i][j] == num and (i, j) != pos:
                return False

    return True

def main():
    # The Sudoku grid parsed from the input. 0 represents an empty cell.
    puzzle = [
        [5, 0, 0, 0, 0, 8, 0, 4, 9],
        [0, 0, 0, 5, 0, 0, 0, 0, 3],
        [0, 0, 6, 0, 7, 3, 0, 0, 1],
        [1, 0, 5, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 8, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 8],
        [7, 0, 0, 0, 0, 4, 1, 0, 5],
        [0, 0, 3, 0, 0, 0, 0, 2, 0],
        [4, 0, 9, 0, 5, 0, 0, 0, 3]
    ]

    # Solve the puzzle
    if solve_sudoku(puzzle):
        # Extract the top horizontal line from the solved puzzle
        top_line = puzzle[0]
        # Print the numbers separated by spaces
        print(' '.join(map(str, top_line)))
    else:
        print("No solution exists for the given Sudoku.")

if __name__ == "__main__":
    main()
