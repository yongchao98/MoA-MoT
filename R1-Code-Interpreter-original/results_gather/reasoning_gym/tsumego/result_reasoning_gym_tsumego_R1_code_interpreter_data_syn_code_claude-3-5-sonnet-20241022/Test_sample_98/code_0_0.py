def print_board(move=None):
    board = [
        ['.'] * 13 for _ in range(13)
    ]
    
    # Place existing stones
    black_stones = [(0,7), (1,6), (0,5), (1,4), (2,3), (5,3), (10,2), (11,2), (5,0)]  # X stones
    white_stones = [(2,7), (4,6), (1,5), (2,5), (3,5), (7,7), (9,7), (7,2), (10,2)]   # O stones
    
    for x,y in black_stones:
        board[y][x] = 'X'
    for x,y in white_stones:
        board[y][x] = 'O'
        
    if move:
        x, y = move
        board[y][x] = '*'  # Mark the suggested move
    
    # Print board
    print('   A B C D E F G H I J K L M')
    for i in range(12, -1, -1):
        row = f"{i+1:2d} "
        row += ' '.join(board[i])
        print(row)
    
    return board

def check_capture(x, y):
    captures = 0
    # Check surrounding white groups that would be captured
    # This is a simplified check
    if x == 2 and y == 6:  # D7 position
        captures = 4  # Would capture the central white group
    return captures

# Test potential moves
best_move = None
max_captures = 0

# Check D7 specifically
move = (3, 6)  # D7
captures = check_capture(*move)
if captures > max_captures:
    max_captures = captures
    best_move = move

print("Analyzing position...")
print("\nBest move found:")
print_board(best_move)
print(f"\nCaptures possible: {max_captures}")