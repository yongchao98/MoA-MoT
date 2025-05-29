# Define the initial positions of the knights
white_knights = {'C1': 'w', 'D1': 'w'}
black_knights = {'C3': 'B', 'B1': 'B'}
empty_squares = {'A1', 'A2', 'B2', 'D2'}

# Define the target positions for the knights
target_white_positions = {'C3', 'B1'}
target_black_positions = {'C1', 'D1'}

# Define the possible moves for a knight from a given position
def knight_moves(position):
    col, row = position
    moves = []
    col_idx = ord(col) - ord('A')
    row_idx = int(row) - 1
    potential_moves = [
        (col_idx + 2, row_idx + 1), (col_idx + 2, row_idx - 1),
        (col_idx - 2, row_idx + 1), (col_idx - 2, row_idx - 1),
        (col_idx + 1, row_idx + 2), (col_idx + 1, row_idx - 2),
        (col_idx - 1, row_idx + 2), (col_idx - 1, row_idx - 2)
    ]
    for c, r in potential_moves:
        if 0 <= c < 4 and 0 <= r < 3:
            moves.append((chr(c + ord('A')), str(r + 1)))
    return moves

# Check if a move is valid
def is_valid_move(start, end, occupied):
    return end in knight_moves(start) and end not in occupied

# Simulate the moves
def simulate_moves():
    moves = []
    occupied = set(white_knights.keys()).union(black_knights.keys())
    occupied = occupied.union(empty_squares)
    
    # Try to find a sequence of moves
    def backtrack(w_positions, b_positions, move_list):
        if set(w_positions.keys()) == target_white_positions and set(b_positions.keys()) == target_black_positions:
            return move_list
        
        # White's turn
        for w_pos in list(w_positions.keys()):
            for move in knight_moves(w_pos):
                if is_valid_move(w_pos, move, occupied):
                    # Make the move
                    w_positions[move] = w_positions.pop(w_pos)
                    occupied.remove(w_pos)
                    occupied.add(move)
                    move_list.append(f"w,{w_pos},{move}")
                    
                    # Black's turn
                    for b_pos in list(b_positions.keys()):
                        for b_move in knight_moves(b_pos):
                            if is_valid_move(b_pos, b_move, occupied):
                                # Make the move
                                b_positions[b_move] = b_positions.pop(b_pos)
                                occupied.remove(b_pos)
                                occupied.add(b_move)
                                move_list.append(f"B,{b_pos},{b_move}")
                                
                                # Recurse
                                result = backtrack(w_positions, b_positions, move_list)
                                if result:
                                    return result
                                
                                # Undo black's move
                                b_positions[b_pos] = b_positions.pop(b_move)
                                occupied.remove(b_move)
                                occupied.add(b_pos)
                                move_list.pop()
                    
                    # Undo white's move
                    w_positions[w_pos] = w_positions.pop(move)
                    occupied.remove(move)
                    occupied.add(w_pos)
                    move_list.pop()
        
        return None
    
    result = backtrack(white_knights.copy(), black_knights.copy(), moves)
    return result

# Get the result
result = simulate_moves()
if result:
    print(result)
else:
    print("No")