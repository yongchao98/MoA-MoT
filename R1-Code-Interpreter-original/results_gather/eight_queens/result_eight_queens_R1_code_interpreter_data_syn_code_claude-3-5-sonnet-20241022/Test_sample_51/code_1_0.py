def print_board(positions, forbidden):
    board = [['0' for _ in range(8)] for _ in range(8)]
    
    # Place queens
    for pos in positions:
        row, col = map(int, pos.split())
        board[row][col] = 'Q'
    
    # Mark forbidden positions
    for row, col in forbidden:
        board[row][col] = 'X'
    
    # Print board and verify
    for row in board:
        print(' '.join(row))
    
    # Verify no queens in forbidden positions
    queens = [(int(pos.split()[0]), int(pos.split()[1])) for pos in positions]
    forbidden = [(4,4), (6,6), (7,4)]
    
    for queen in queens:
        if queen in forbidden:
            return False
            
    return True

positions = ['0 0', '1 6', '2 4', '3 7', '4 1', '5 3', '6 5', '7 2']
forbidden = [(4,4), (6,6), (7,4)]
print("Board configuration:")
is_valid = print_board(positions, forbidden)
print("\nValid solution:", is_valid)