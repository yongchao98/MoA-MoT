def verify_final_moves():
    # Start from the state we know is valid
    board = {
        (4,2): 'B', (2,1): 'w', (3,3): 'w', (1,1): 'B'
    }
    
    # Remaining moves to verify
    remaining_moves = ['w3341', 'B1132', 'w4122']
    
    for move in remaining_moves:
        color = move[0]
        from_x, from_y = int(move[1]), int(move[2])
        to_x, to_y = int(move[3]), int(move[4])
        
        # Make the move
        if (from_x, from_y) in board:
            board[(to_x, to_y)] = board[(from_x, from_y)]
            del board[(from_x, from_y)]
            print(f"After move {move}:", board)
    
    # Check final positions
    target = {
        (2,1): 'w', (2,2): 'w',  # B1,B2
        (3,2): 'B', (4,2): 'B'   # C2,D2
    }
    
    print("Final board:", board)
    print("Target board:", target)
    return board == target

print("Final moves valid:", verify_final_moves())