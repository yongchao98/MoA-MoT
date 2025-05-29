def is_valid_move(start, end):
    # Calculate the difference in coordinates
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights():
    # Initial positions
    white_knights = {'A1', 'D2'}
    black_knights = {'A3', 'B1'}
    empty_squares = {'C1', 'B3', 'C2', 'C3', 'D1', 'D3'}
    
    # Target positions
    target_white = {'A3', 'B1'}
    target_black = {'A1', 'D2'}
    
    # Moves list
    moves = []
    
    # Function to perform a move
    def move_knight(knight_set, from_pos, to_pos):
        knight_set.remove(from_pos)
        knight_set.add(to_pos)
        empty_squares.remove(to_pos)
        empty_squares.add(from_pos)
        moves.append(f"{'w' if from_pos in white_knights else 'B'},{from_pos},{to_pos}")
    
    # Try to find a sequence of moves
    def find_moves(turn):
        # Check if the current configuration is the target configuration
        if white_knights == target_white and black_knights == target_black:
            return True
        
        # Determine whose turn it is
        current_knights = white_knights if turn == 'w' else black_knights
        
        # Try all possible moves for the current player's knights
        for knight in list(current_knights):
            for es in list(empty_squares):
                if is_valid_move(knight, es):
                    move_knight(current_knights, knight, es)
                    if find_moves('B' if turn == 'w' else 'w'):
                        return True
                    # Backtrack
                    move_knight(current_knights, es, knight)
        
        return False
    
    if find_moves('w'):
        return moves
    else:
        return "No"

# Execute the function and print the result
print(swap_knights())