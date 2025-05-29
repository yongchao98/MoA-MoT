def get_knight_moves(x, y):
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(nx, ny) for nx, ny in moves if 1 <= nx <= 4 and 1 <= ny <= 3]

def verify_impossible():
    # Convert board coordinates to numerical (A=1, B=2, C=3, D=4)
    board = {
        (3,3): 'B', (4,3): 'w',  # C3, D3
        (3,1): 'B', (4,1): 'w',  # C1, D1
        (1,2): '.', (2,2): '.'   # A2, B2
    }
    
    # Check reachability from initial positions
    white_positions = [(4,3), (4,1)]
    black_positions = [(3,3), (3,1)]
    empty_squares = [(1,2), (2,2)]
    
    # Print all possible first moves for white
    print("Possible first moves for white knights:")
    for pos in white_positions:
        moves = get_knight_moves(pos[0], pos[1])
        valid_moves = [m for m in moves if m in empty_squares]
        print(f"From {pos}: {valid_moves}")
    
    # Print all possible first moves for black after any white move
    print("\nPossible first moves for black knights:")
    for pos in black_positions:
        moves = get_knight_moves(pos[0], pos[1])
        valid_moves = [m for m in moves if m in empty_squares]
        print(f"From {pos}: {valid_moves}")

verify_impossible()