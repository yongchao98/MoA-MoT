def is_valid_move(start, end):
    # Calculate the difference in coordinates
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = {'A3', 'D3'}
    black_knights = {'B3', 'A1'}
    empty_squares = {'C3', 'D2', 'B1', 'C1'}
    
    # Target positions
    white_target = {'B3', 'A1'}
    black_target = {'A3', 'D3'}
    
    # Moves list
    moves = []
    
    # Function to perform the move
    def move_knight(color, from_pos, to_pos):
        if color == 'w':
            white_knights.remove(from_pos)
            white_knights.add(to_pos)
        else:
            black_knights.remove(from_pos)
            black_knights.add(to_pos)
        empty_squares.remove(to_pos)
        empty_squares.add(from_pos)
        moves.append(f"{color},{from_pos},{to_pos}")
    
    # Backtracking function to find the solution
    def backtrack(turn):
        if white_knights == white_target and black_knights == black_target:
            return True
        
        if turn == 'w':
            for wk in list(white_knights):
                for es in list(empty_squares):
                    if is_valid_move(wk, es):
                        move_knight('w', wk, es)
                        if backtrack('b'):
                            return True
                        # Undo move
                        move_knight('w', es, wk)
                        moves.pop()
        else:
            for bk in list(black_knights):
                for es in list(empty_squares):
                    if is_valid_move(bk, es):
                        move_knight('b', bk, es)
                        if backtrack('w'):
                            return True
                        # Undo move
                        move_knight('b', es, bk)
                        moves.pop()
        
        return False
    
    if backtrack('w'):
        return moves
    else:
        return "No"

# Find and print the moves
print(find_moves())