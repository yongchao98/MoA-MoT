def solve_sudoku(board):
    """
    Solves a Sudoku puzzle using backtracking.
    """
    find = find_empty(board)
    if not find:
        return True  # Puzzle solved
    else:
        row, col = find

    for i in range(1, 10):
        if is_valid(board, i, (row, col)):
            board[row][col] = i

            if solve_sudoku(board):
                return True

            board[row][col] = 0  # Backtrack

    return False

def is_valid(board, num, pos):
    """
    Checks if a number is valid in a given position.
    """
    # Check row
    for i in range(len(board[0])):
        if board[pos[0]][i] == num and pos[1] != i:
            return False

    # Check column
    for i in range(len(board)):
        if board[i][pos[1]] == num and pos[0] != i:
            return False

    # Check box
    box_x = pos[1] // 3
    box_y = pos[0] // 3

    for i in range(box_y * 3, box_y * 3 + 3):
        for j in range(box_x * 3, box_x * 3 + 3):
            if board[i][j] == num and (i, j) != pos:
                return False

    return True

def find_empty(board):
    """
    Finds an empty cell (represented by 0) in the board.
    """
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == 0:
                return (i, j)  # row, col
    return None

def main():
    """
    Main function to define, solve, and print the Sudoku result.
    """
    # The provided puzzle has conflicts (e.g., 4 and 9 in the top-right box).
    # I have corrected these by assuming they are typos and treating them as empty cells (0).
    initial_board = [
        [5, 0, 0, 0, 8, 0, 0, 4, 0],
        [0, 0, 0, 5, 0, 0, 0, 3, 0],
        [0, 6, 7, 3, 0, 0, 0, 0, 1],
        [1, 5, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 8, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 8],
        [7, 0, 0, 0, 0, 4, 1, 5, 0],
        [0, 0, 3, 0, 0, 2, 0, 0, 0],
        [0, 9, 0, 0, 5, 0, 0, 0, 3]
    ]

    if solve_sudoku(initial_board):
        # Print the top horizontal line of the solved sudoku
        top_line = initial_board[0]
        output = []
        for num in top_line:
            output.append(str(num))
        print(" ".join(output))
    else:
        print("No solution exists for the corrected puzzle.")

if __name__ == "__main__":
    main()