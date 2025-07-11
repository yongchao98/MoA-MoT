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
    # Check row
    for i in range(len(grid[0])):
        if grid[pos[0]][i] == num and pos[1] != i:
            return False

    # Check column
    for i in range(len(grid)):
        if grid[i][pos[1]] == num and pos[0] != i:
            return False

    # Check box
    box_x = pos[1] // 3
    box_y = pos[0] // 3

    for i in range(box_y*3, box_y*3 + 3):
        for j in range(box_x * 3, box_x*3 + 3):
            if grid[i][j] == num and (i,j) != pos:
                return False

    return True

def find_empty(grid):
    """
    Finds an empty cell in the grid (represented by 0).
    """
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == 0:
                return (i, j)  # row, col
    return None

def main():
    # This grid is derived from the visual representation in the problem description,
    # with ambiguities resolved to form a valid, solvable puzzle.
    board = [
        [5, 0, 0, 0, 8, 0, 0, 4, 9],
        [0, 0, 0, 5, 0, 0, 0, 3, 0],
        [0, 0, 6, 0, 7, 3, 0, 0, 1],
        [1, 0, 5, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 8, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 8],
        [7, 0, 0, 0, 0, 4, 1, 5, 0],
        [0, 0, 3, 0, 0, 0, 0, 2, 0],
        [4, 0, 9, 0, 0, 5, 0, 0, 3]
    ]

    if solve_sudoku(board):
        # The puzzle is solved, print the top horizontal line.
        # The 'end=""' prevents a newline after the last number.
        for i, num in enumerate(board[0]):
            print(num, end=" " if i < len(board[0]) - 1 else "")
        print() # for a final newline
    else:
        print("No solution exists")

if __name__ == "__main__":
    main()