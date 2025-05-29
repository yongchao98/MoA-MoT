def is_safe(board, row, col):
    # Check row
    for j in range(8):
        if board[row][j] == 1:
            return False
    
    # Check column
    for i in range(8):
        if board[i][col] == 1:
            return False
    
    # Check diagonals
    for i in range(8):
        for j in range(8):
            if board[i][j] == 1:
                if abs(row - i) == abs(col - j):
                    return False
    
    return True

def print_solution(board):
    positions = []
    for i in range(8):
        for j in range(8):
            if board[i][j] == 1:
                positions.append(f"{i} {j}")
    print(",".join(positions))

def solve_queens():
    # Initialize board with constraints
    board = [[0 for x in range(8)] for y in range(8)]
    # Set initial constraints
    board[0][5] = 1  # Pre-placed queen
    forbidden = [(0,6), (3,0), (5,0), (7,3)]
    for x, y in forbidden:
        board[x][y] = 'X'
    
    def solve_util(queens_placed):
        if queens_placed == 8:
            return True
        
        # Try placing queen in each position
        for col in range(8):
            if col == 5:  # Skip column with pre-placed queen
                continue
            for row in range(8):
                if board[row][col] == 0 and is_safe(board, row, col):
                    # Try placing queen
                    board[row][col] = 1
                    if solve_util(queens_placed + 1):
                        return True
                    board[row][col] = 0
        
        return False
    
    if solve_util(1):  # Start with 1 queen already placed
        print_solution(board)
    else:
        print("No solution exists")

solve_queens()