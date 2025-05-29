def is_safe(row, col, queens, blocked):
    # Check if the column is safe
    for r, c in queens:
        if c == col or abs(r - row) == abs(c - col):
            return False
    return True

def solve_n_queens(n, blocked):
    def backtrack(row, queens):
        if row == n:
            return queens
        for col in range(n):
            if (row, col) not in blocked and is_safe(row, col, queens, blocked):
                queens.append((row, col))
                result = backtrack(row + 1, queens)
                if result:
                    return result
                queens.pop()
        return None

    blocked_positions = set()
    for r in range(n):
        for c in range(n):
            if blocked[r][c] == 'X':
                blocked_positions.add((r, c))

    return backtrack(0, [])

# Define the board with blocked positions
blocked_board = [
    "00000000",
    "00000000",
    "00000000",
    "00000000",
    "00000000",
    "XXX00000",
    "00000000",
    "00000000"
]

# Find the solution
solution = solve_n_queens(8, blocked_board)
print(solution)