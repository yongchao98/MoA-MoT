def is_safe(board, row, col, n):
    # Check row
    for j in range(n):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(n):
        if board[i][col] == 1:
            return False
    
    # Check diagonals
    for i in range(n):
        for j in range(n):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def solve_n_queens(n):
    # Initialize board with pre-placed queens
    board = [
        [0,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0]
    ]
    
    queens_positions = [(1,0), (3,5), (6,1)]  # Pre-placed queens
    remaining_rows = [i for i in range(n) if i not in [1,3,6]]
    
    def solve(queens_left):
        if queens_left == 0:
            return True
            
        for row in remaining_rows:
            for col in range(n):
                if is_safe(board, row, col, n):
                    board[row][col] = 1
                    queens_positions.append((row, col))
                    if solve(queens_left - 1):
                        return True
                    board[row][col] = 0
                    queens_positions.pop()
        return False
    
    if solve(5):  # We need to place 5 more queens
        # Sort positions and format output
        queens_positions.sort()
        result = ", ".join(f"{pos[0]} {pos[1]}" for pos in queens_positions)
        print(f"<<<{result}>>>")
    else:
        print("No solution exists!")

solve_n_queens(8)