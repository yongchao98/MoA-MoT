import collections

def get_rook_attacks(r, c, occupied):
    """ Get all squares attacked by a rook. """
    attacks = set()
    # Horizontal right
    for i in range(c + 1, 8):
        attacks.add((r, i))
        if (r, i) in occupied: break
    # Horizontal left
    for i in range(c - 1, -1, -1):
        attacks.add((r, i))
        if (r, i) in occupied: break
    # Vertical down
    for i in range(r + 1, 8):
        attacks.add((i, c))
        if (i, c) in occupied: break
    # Vertical up
    for i in range(r - 1, -1, -1):
        attacks.add((i, c))
        if (i, c) in occupied: break
    return attacks

def get_bishop_attacks(r, c, occupied):
    """ Get all squares attacked by a bishop. """
    attacks = set()
    # Down-right
    for i in range(1, 8):
        pos = (r + i, c + i)
        if 0 <= pos[0] < 8 and 0 <= pos[1] < 8:
            attacks.add(pos)
            if pos in occupied: break
        else: break
    # Down-left
    for i in range(1, 8):
        pos = (r + i, c - i)
        if 0 <= pos[0] < 8 and 0 <= pos[1] < 8:
            attacks.add(pos)
            if pos in occupied: break
        else: break
    # Up-right
    for i in range(1, 8):
        pos = (r - i, c + i)
        if 0 <= pos[0] < 8 and 0 <= pos[1] < 8:
            attacks.add(pos)
            if pos in occupied: break
        else: break
    # Up-left
    for i in range(1, 8):
        pos = (r - i, c - i)
        if 0 <= pos[0] < 8 and 0 <= pos[1] < 8:
            attacks.add(pos)
            if pos in occupied: break
        else: break
    return attacks

def get_knight_moves(r, c):
    """ Get all squares attacked by a knight. """
    attacks = set()
    moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
             (1, -2), (1, 2), (2, -1), (2, 1)]
    for dr, dc in moves:
        if 0 <= r + dr < 8 and 0 <= c + dc < 8:
            attacks.add((r + dr, c + dc))
    return attacks

def main():
    """
    Main function to verify the chess puzzle solution.
    """
    # Board representation: (row, col) from (0,0) [a8] to (7,7) [h1]
    # Files: a=0, b=1, ..., h=7
    # Ranks: 8=0, 7=1, ..., 1=7
    
    # Proposed pieces and their values
    pieces = {
        'Q': {'pos': (6, 3), 'val': 9},  # Queen on d2
        'N': {'pos': (5, 5), 'val': 3},  # Knight on f3
        'B': {'pos': (7, 2), 'val': 3},  # Bishop on c1
    }
    
    black_king_pos = (7, 7) # King on h1
    
    occupied_squares = {p['pos'] for p in pieces.values()}
    
    total_attacked_squares = set()

    for piece_type, piece_info in pieces.items():
        r, c = piece_info['pos']
        if piece_type == 'Q':
            total_attacked_squares.update(get_rook_attacks(r, c, occupied_squares))
            total_attacked_squares.update(get_bishop_attacks(r, c, occupied_squares))
        elif piece_type == 'N':
            total_attacked_squares.update(get_knight_moves(r, c))
        elif piece_type == 'B':
            total_attacked_squares.update(get_bishop_attacks(r, c, occupied_squares))
            
    # Verification
    king_is_safe = black_king_pos not in total_attacked_squares
    
    king_adjacent_squares = {
        (black_king_pos[0] - 1, black_king_pos[1] - 1), # g2
        (black_king_pos[0] - 1, black_king_pos[1]),     # h2
        (black_king_pos[0],     black_king_pos[1] - 1)  # g1
    }
    
    stalemate_condition_met = king_adjacent_squares.issubset(total_attacked_squares)
    
    num_unattacked = 64 - len(total_attacked_squares)
    only_one_unattacked = num_unattacked == 1
    
    print("Verifying the proposed solution...")
    print(f"White pieces: Queen at d2, Knight at f3, Bishop at c1")
    print(f"Black king at h1 is stalemated.")
    print("-" * 20)
    print(f"Is king's square (h1) unattacked? {'Yes' if king_is_safe else 'No'}")
    print(f"Are all adjacent squares attacked? {'Yes' if stalemate_condition_met else 'No'}")
    print(f"Total squares attacked: {len(total_attacked_squares)}")
    print(f"Total squares unattacked: {num_unattacked}")
    print(f"Is there exactly one unattacked square on the board? {'Yes' if only_one_unattacked else 'No'}")
    print("-" * 20)

    if king_is_safe and stalemate_condition_met and only_one_unattacked:
        print("Verification successful!")
        queen_val = pieces['Q']['val']
        knight_val = pieces['N']['val']
        bishop_val = pieces['B']['val']
        total_val = queen_val + knight_val + bishop_val

        # Final print statement showing the breakdown as requested
        print("The calculation for the smallest number of points is:")
        print("Queen (", queen_val, ") + Knight (", knight_val, ") + Bishop (", bishop_val, ") = ", total_val, sep="")

    else:
        print("The proposed solution is incorrect.")

if __name__ == "__main__":
    main()
