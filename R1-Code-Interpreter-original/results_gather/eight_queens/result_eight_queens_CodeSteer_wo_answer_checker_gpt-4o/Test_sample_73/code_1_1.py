def solve_queens(board):
    n = 8
    queens = []
    occupied_cols = set()
    occupied_major_diags = set()
    occupied_minor_diags = set()

    # Mark already placed queens and blocked positions
    for r in range(n):
        for c in range(n):
            if board[r][c] == '1':
                queens.append((r, c))
                occupied_cols.add(c)
                occupied_major_diags.add(r - c)
                occupied_minor_diags.add(r + c)

    # Backtracking function to place queens
    def place_queens(row):
        if row == n:
            return True  # All queens are placed
        for col in range(n):
            if (col not in occupied_cols and
                (row - col) not in occupied_major_diags and
                (row + col) not in occupied_minor_diags and
                board[row][col] != 'X'):
                
                # Place queen
                queens.append((row, col))
                occupied_cols.add(col)
                occupied_major_diags.add(row - col)
                occupied_minor_diags.add(row + col)

                # Recurse to place the next queen
                if place_queens(row + 1):
                    return True

                # Backtrack
                queens.pop()
                occupied_cols.remove(col)
                occupied_major_diags.remove(row - col)
                occupied_minor_diags.remove(row + col)
        return False

    # Start placing queens from the first row
    place_queens(0)

    # Format the output
    result = ', '.join(f'{r} {c}' for r, c in sorted(queens))
    return f'<<<{result}>>>'

# Define the board
board = [
    "00000000",
    "00010000",
    "00000000",
    "00000000",
    "0000000X",
    "00X00000",
    "00000000",
    "X0000000"
]

# Solve and print the result
print(solve_queens(board))