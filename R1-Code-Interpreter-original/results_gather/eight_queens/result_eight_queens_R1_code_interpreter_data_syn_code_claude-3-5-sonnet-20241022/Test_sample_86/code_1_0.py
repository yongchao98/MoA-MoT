def is_safe(board, row, col, n):
    # Check row on left side and right side
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column on upper and lower side
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonals in all four directions
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_queens(board, col, n, blocked_pos, placed_queen):
    # Base case: If all queens are placed, return True
    if col >= n:
        return True
    
    # Skip the column where queen is already placed
    if col == placed_queen[1]:
        return solve_queens(board, col + 1, n, blocked_pos, placed_queen)
    
    # Consider this column and try placing queens in all rows one by one
    for row in range(n):
        # Check if we can place a queen here
        if (row, col) != blocked_pos and board[row][col] != 1:
            # Check if position is safe considering all placed queens
            temp = board[row][col]
            board[row][col] = 1
            if is_safe(board, row, col, n):
                # Recur to place rest of the queens
                if solve_queens(board, col + 1, n, blocked_pos, placed_queen):
                    return True
            # Backtrack
            board[row][col] = temp
    
    return False

def find_queens_positions():
    n = 8
    # Initialize board with given state
    board = [[0]*n for _ in range(n)]
    placed_queen = (2, 6)  # Pre-placed queen position
    board[placed_queen[0]][placed_queen[1]] = 1
    blocked_pos = (6, 3)
    
    # Start from column 0
    if solve_queens(board, 0, n, blocked_pos, placed_queen):
        # Collect positions of queens
        positions = []
        for i in range(n):
            for j in range(n):
                if board[i][j] == 1:
                    positions.append(f"{i} {j}")
        print("<<<" + ", ".join(sorted(positions)) + ">>>")
    else:
        print("No solution exists")

find_queens_positions()