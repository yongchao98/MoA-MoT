# Function to check if a move at (x, y) can capture any stones
def can_capture(x, y):
    # Check if the move is within bounds and the position is empty
    if 0 <= x < 13 and 0 <= y < 13 and board[y][x] == '.':
        # Simulate the move
        board[y][x] = 'X'
        # Check surrounding positions for potential captures
        # This is a simplified check for demonstration purposes
        # In a real game, we would need to check for liberties and full capture
        if (y > 0 and board[y-1][x] == 'O') or (y < 12 and board[y+1][x] == 'O') or \
           (x > 0 and board[y][x-1] == 'O') or (x < 12 and board[y][x+1] == 'O'):
            return True
        # Revert the move
        board[y][x] = '.'
    return False

# Check the potential move at B3 (x=1, y=10)
move_x, move_y = 1, 10
capture_possible = can_capture(move_x, move_y)
print(capture_possible)