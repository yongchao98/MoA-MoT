def is_safe(board, row, col, n):
    # Check this row
    for j in range(n):
        if board[row][j] == 1:
            return False
            
    # Check this column
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

def solve_n_queens(n=8):
    # Initialize board with the given queen at (3,4)
    board = [[0 for x in range(n)] for y in range(n)]
    board[3][4] = 1  # Given queen
    queens_placed = [(3,4)]  # Keep track of placed queens
    
    def solve_util(queens_to_place):
        if queens_to_place == 0:
            return True
            
        for row in range(n):
            for col in range(n):
                # Skip if this position already has a queen
                if (row, col) in queens_placed:
                    continue
                    
                # Try placing a queen if position is empty and safe
                if board[row][col] == 0:
                    # Temporarily place queen and check if it's safe
                    board[row][col] = 1
                    is_position_safe = True
                    
                    # Check against all placed queens
                    for queen_row, queen_col in queens_placed:
                        if (queen_row == row or queen_col == col or 
                            abs(queen_row - row) == abs(queen_col - col)):
                            is_position_safe = False
                            break
                    
                    if is_position_safe:
                        queens_placed.append((row, col))
                        if solve_util(queens_to_place - 1):
                            return True
                        queens_placed.pop()
                    board[row][col] = 0
        return False
    
    if solve_util(7):  # We need to place 7 more queens
        # Sort positions for consistent output
        queens_placed.sort()
        positions = [f"{row} {col}" for row, col in queens_placed]
        print(",".join(positions))
    else:
        print("No solution exists")

solve_n_queens()