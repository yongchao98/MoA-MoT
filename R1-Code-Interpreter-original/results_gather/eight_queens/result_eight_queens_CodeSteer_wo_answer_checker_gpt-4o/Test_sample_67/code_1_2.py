def is_safe(row, col, cols, major_diags, minor_diags):
    return col not in cols and (row - col) not in major_diags and (row + col) not in minor_diags

def solve_n_queens(board, row, cols, major_diags, minor_diags, queens_positions):
    if row >= len(board):
        return True

    for col in range(len(board)):
        if board[row][col] == 0 and is_safe(row, col, cols, major_diags, minor_diags):
            # Place the queen
            board[row][col] = 1
            cols.add(col)
            major_diags.add(row - col)
            minor_diags.add(row + col)
            queens_positions.append((row, col))

            # Move to the next row
            if solve_n_queens(board, row + 1, cols, major_diags, minor_diags, queens_positions):
                return True

            # Backtrack
            board[row][col] = 0
            cols.remove(col)
            major_diags.remove(row - col)
            minor_diags.remove(row + col)
            queens_positions.pop()

    return False

def main():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0]
    ]

    # Initial queen position
    queens_positions = [(7, 1)]

    # Track columns and diagonals under attack
    cols = {1}
    major_diags = {7 - 1}
    minor_diags = {7 + 1}

    # Start solving from the first row
    if solve_n_queens(board, 0, cols, major_diags, minor_diags, queens_positions):
        result = ', '.join(f"{r} {c}" for r, c in queens_positions)
        print(f"<<<{result}>>>")
    else:
        print("No solution found")

main()