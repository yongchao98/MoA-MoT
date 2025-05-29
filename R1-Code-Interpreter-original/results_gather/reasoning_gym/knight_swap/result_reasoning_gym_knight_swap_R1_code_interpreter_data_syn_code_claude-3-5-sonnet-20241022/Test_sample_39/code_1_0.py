def is_valid_knight_move(from_pos, to_pos):
    dx = abs(from_pos[0] - to_pos[0])
    dy = abs(from_pos[1] - to_pos[1])
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def verify_solution(moves):
    # Initial board
    board = {
        (2,1): 'B', (2,2): 'B',  # B1,B2
        (3,2): 'w', (4,2): 'w'   # C2,D2
    }
    
    print("Initial board:", board)
    
    for move in moves:
        color = move[0]
        from_x, from_y = int(move[1]), int(move[2])
        to_x, to_y = int(move[3]), int(move[4])
        
        # Check if piece exists
        if (from_x, from_y) not in board or board[(from_x, from_y)] != color:
            print(f"Invalid move {move}: No {color} piece at {from_x},{from_y}")
            return False
            
        # Check if destination is empty
        if (to_x, to_y) in board:
            print(f"Invalid move {move}: Destination {to_x},{to_y} is occupied")
            return False
            
        # Check if move is valid knight move
        if not is_valid_knight_move((from_x, from_y), (to_x, to_y)):
            print(f"Invalid move {move}: Not a valid knight move")
            return False
            
        # Make the move
        board[(to_x, to_y)] = board[(from_x, from_y)]
        del board[(from_x, from_y)]
        print(f"After move {move}:", board)
    
    # Check final positions
    target = {
        (2,1): 'w', (2,2): 'w',  # B1,B2
        (3,2): 'B', (4,2): 'B'   # C2,D2
    }
    
    if board == target:
        return True
    else:
        print("Final positions don't match target")
        return False

moves = ['B2133', 'w3213', 'B2243', 'w4223', 'B3321', 'w2331', 'B2142', 'w3112', 
         'B4331', 'w1321', 'B3123', 'w1233', 'B2311', 'w3341', 'B1132', 'w4122']

print("Solution valid:", verify_solution(moves))