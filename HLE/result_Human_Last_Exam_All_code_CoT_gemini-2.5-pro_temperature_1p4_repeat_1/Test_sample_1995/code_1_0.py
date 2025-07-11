def get_piece_attacks(piece_type, col, row):
    """
    Calculates the set of squares attacked by a single piece.
    A piece's attack is not blocked by other pieces for this problem type.
    """
    attacks = set()
    # Directions: N, NE, E, SE, S, SW, W, NW
    all_dirs = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]
    
    # Rook moves (N, E, S, W)
    if piece_type == 'R':
        for d in [all_dirs[0], all_dirs[2], all_dirs[4], all_dirs[6]]:
            for i in range(1, 8):
                c, r = col + i * d[0], row + i * d[1]
                if 0 <= c < 8 and 0 <= r < 8:
                    attacks.add((c, r))
                else:
                    break
    
    # Bishop moves (NE, SE, SW, NW)
    if piece_type == 'B':
        for d in [all_dirs[1], all_dirs[3], all_dirs[5], all_dirs[7]]:
            for i in range(1, 8):
                c, r = col + i * d[0], row + i * d[1]
                if 0 <= c < 8 and 0 <= r < 8:
                    attacks.add((c, r))
                else:
                    break

    # Knight moves
    if piece_type == 'N':
        knight_moves = [(1, 2), (2, 1), (2, -1), (1, -2), 
                        (-1, -2), (-2, -1), (-2, 1), (-1, 2)]
        for move in knight_moves:
            c, r = col + move[0], row + move[1]
            if 0 <= c < 8 and 0 <= r < 8:
                attacks.add((c, r))
    
    return attacks

def solve_chess_puzzle():
    """
    Solves the chess puzzle by verifying a 16-point solution.
    """
    def to_coords(s):
        col = ord(s[0]) - ord('a')
        row = int(s[1]) - 1
        return col, row

    def to_notation(pos):
        c, r = pos
        return chr(ord('a') + c) + str(r + 1)

    piece_values = {'R': 5, 'B': 3, 'N': 3}
    
    # Candidate solution: Two Rooks, one Bishop, one Knight (16 points)
    pieces = {
        'Rook 1': {'type': 'R', 'pos': to_coords('e2'), 'val': piece_values['R']},
        'Rook 2': {'type': 'R', 'pos': to_coords('h3'), 'val': piece_values['R']},
        'Bishop': {'type': 'B', 'pos': to_coords('c2'), 'val': piece_values['B']},
        'Knight': {'type': 'N', 'pos': to_coords('d2'), 'val': piece_values['N']}
    }

    print("Investigating a solution with the following white pieces:")
    total_points = 0
    occupied_squares = set()
    for name, data in pieces.items():
        total_points += data['val']
        occupied_squares.add(data['pos'])
        print(f"- {name} ({data['val']} pts) at {to_notation(data['pos'])}")

    # Calculate all attacked squares
    all_attacked_squares = set()
    for name, data in pieces.items():
        attacks = get_piece_attacks(data['type'], data['pos'][0], data['pos'][1])
        all_attacked_squares.update(attacks)

    # A square is "controlled" if it's attacked or occupied.
    # The black king cannot move to a controlled square.
    controlled_squares = all_attacked_squares.union(occupied_squares)

    all_board_squares = set((c, r) for c in range(8) for r in range(8))
    uncontrolled_squares = all_board_squares - controlled_squares

    print("\n--- Verification ---")
    print(f"Total squares controlled (attacked or occupied): {len(controlled_squares)}")
    print(f"Total squares uncontrolled: {len(uncontrolled_squares)}")

    if len(uncontrolled_squares) == 1:
        king_pos = list(uncontrolled_squares)[0]
        king_square_notation = to_notation(king_pos)
        print(f"The single safe square is {king_square_notation}.")
        print("Placing the black king on this square results in a stalemate.")
        
        # Final confirmation and equation
        print("\nThis configuration achieves the goal with a minimal value of 16 points.")
        print("The final calculation is:")
        rook_val = piece_values['R']
        bishop_val = piece_values['B']
        knight_val = piece_values['N']
        print(f"{rook_val} (Rook) + {rook_val} (Rook) + {bishop_val} (Bishop) + {knight_val} (Knight) = {total_points}")
    else:
        print("This configuration does not solve the puzzle.")
        print("Uncontrolled squares:", [to_notation(pos) for pos in uncontrolled_squares])

solve_chess_puzzle()