def is_safe(row, col, cols, diag1, diag2):
    return col not in cols and (row - col) not in diag1 and (row + col) not in diag2

def solve_n_queens(n, board):
    def backtrack(row, cols, diag1, diag2, positions):
        if row == n:
            return positions
        for col in range(n):
            if board[row][col] == 'X':
                continue
            if is_safe(row, col, cols, diag1, diag2):
                cols.add(col)
                diag1.add(row - col)
                diag2.add(row + col)
                positions.append((row, col))
                result = backtrack(row + 1, cols, diag1, diag2, positions)
                if result:
                    return result
                positions.pop()
                cols.remove(col)
                diag1.remove(row - col)
                diag2.remove(row + col)
        return None

    cols = set()
    diag1 = set()
    diag2 = set()
    positions = []

    # Mark the initial queen position
    initial_queen = (6, 2)
    cols.add(initial_queen[1])
    diag1.add(initial_queen[0] - initial_queen[1])
    diag2.add(initial_queen[0] + initial_queen[1])
    positions.append(initial_queen)

    return backtrack(0, cols, diag1, diag2, positions)

# Define the board
board = [
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0'],
    ['0', '0', '1', '0', '0', '0', '0', '0'],
    ['0', '0', '0', '0', '0', '0', '0', '0']
]

# Solve the problem
solution = solve_n_queens(8, board)
print(solution)