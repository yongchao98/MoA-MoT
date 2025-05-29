def print_board(positions):
    board = [['0' for _ in range(8)] for _ in range(8)]
    # Mark forbidden positions
    board[2][6] = 'X'
    board[6][2] = 'X'
    board[6][3] = 'X'
    
    # Place queens
    for pos in positions:
        row, col = map(int, pos.split())
        board[row][col] = '1'
    
    # Print board
    print("Board visualization:")
    for row in board:
        print(' '.join(row))
    
    # Verify solution
    queens = [(int(pos.split()[0]), int(pos.split()[1])) for pos in positions]
    
    # Check number of queens
    if len(queens) != 8:
        print(f"Error: Found {len(queens)} queens instead of 8")
        return False
        
    # Check for attacks
    for i, (r1, c1) in enumerate(queens):
        for j, (r2, c2) in enumerate(queens):
            if i != j:
                if (r1 == r2 or  # same row
                    c1 == c2 or  # same column
                    abs(r1 - r2) == abs(c1 - c2)):  # same diagonal
                    print(f"Error: Queens at {r1},{c1} and {r2},{c2} can attack each other")
                    return False
    
    # Check forbidden positions
    forbidden = [(2,6), (6,2), (6,3)]
    for queen in queens:
        if queen in forbidden:
            print(f"Error: Queen placed in forbidden position {queen}")
            return False
            
    print("Solution is valid!")
    return True

positions = "0 1, 1 3, 2 5, 3 7, 4 2, 5 0, 6 6, 7 4"
print_board(positions.split(", "))