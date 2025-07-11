def solve_sudoku(grid):
    """
    Solves a Sudoku puzzle using backtracking.
    """
    find = find_empty(grid)
    if not find:
        return True
    else:
        row, col = find

    for i in range(1, 10):
        if is_valid(grid, i, (row, col)):
            grid[row][col] = i

            if solve_sudoku(grid):
                return True

            grid[row][col] = 0

    return False

def is_valid(grid, num, pos):
    """
    Checks if a number is valid in a given position.
    """
    row, col = pos

    # Check row
    for i in range(len(grid[0])):
        if grid[row][i] == num and col != i:
            return False

    # Check column
    for i in range(len(grid)):
        if grid[i][col] == num and row != i:
            return False

    # Check 3x3 box
    box_x = col // 3
    box_y = row // 3

    for i in range(box_y*3, box_y*3 + 3):
        for j in range(box_x*3, box_x*3 + 3):
            if grid[i][j] == num and (i,j) != pos:
                return False

    return True

def find_empty(grid):
    """
    Finds an empty cell (represented by 0) in the grid.
    """
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == 0:
                return (i, j)  # row, col
    return None

def main():
    """
    Main function to solve the puzzle and print the result.
    """
    # The puzzle is parsed from the visual representation.
    # Note: The original puzzle is invalid with two '5's in the 5th column.
    # The '5' in the last row has been removed to make it solvable.
    grid = [
        [0, 5, 0, 0, 8, 0, 0, 4, 9],
        [0, 0, 0, 0, 5, 0, 0, 3, 0],
        [0, 0, 6, 7, 3, 0, 0, 0, 1],
        [1, 0, 5, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 2, 8, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 8],
        [7, 0, 0, 0, 0, 4, 1, 5, 0],
        [0, 0, 3, 0, 0, 2, 0, 0, 0],
        [4, 0, 9, 0, 0, 0, 0, 0, 3]  # Corrected from original: removed a '5'
    ]

    if solve_sudoku(grid):
        # Print the top horizontal line of the solved sudoku
        solved_first_row = grid[0]
        print(' '.join(map(str, solved_first_row)))
    else:
        print("Could not solve the Sudoku puzzle.")

if __name__ == '__main__':
    main()