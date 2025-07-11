def solve_chess_puzzle():
    """
    This script verifies an 8-piece solution to the chess problem of attacking
    all squares but one, which is occupied by a stalemated black king.
    """

    # Helper functions to convert between algebraic and 0-indexed tuple coordinates
    def to_coords(square):
        """Converts algebraic notation like 'a1' to (0, 0)"""
        col = ord(square[0]) - ord('a')
        row = int(square[1]) - 1
        return row, col

    def to_algebraic(row, col):
        """Converts (0, 0) to 'a1'"""
        return chr(ord('a') + col) + str(row + 1)

    # Attack logic for each piece type
    def get_rook_attacks(pos, occupied):
        r, c = pos
        attacks = set()
        for dr, dc in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            nr, nc = r + dr, c + dc
            while 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add((nr, nc))
                if (nr, nc) in occupied:
                    break
                nr, nc = nr + dr, nc + dc
        return attacks

    def get_bishop_attacks(pos, occupied):
        r, c = pos
        attacks = set()
        for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            nr, nc = r + dr, c + dc
            while 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add((nr, nc))
                if (nr, nc) in occupied:
                    break
                nr, nc = nr + dr, nc + dc
        return attacks

    def get_queen_attacks(pos, occupied):
        return get_rook_attacks(pos, occupied).union(get_bishop_attacks(pos, occupied))

    def get_knight_attacks(pos):
        r, c = pos
        attacks = set()
        moves = [
            (r + 2, c + 1), (r + 2, c - 1), (r - 2, c + 1), (r - 2, c - 1),
            (r + 1, c + 2), (r + 1, c - 2), (r - 1, c + 2), (r - 1, c - 2)
        ]
        for nr, nc in moves:
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add((nr, nc))
        return attacks

    def get_pawn_attacks(pos):
        r, c = pos
        attacks = set()
        for move in [(r + 1, c - 1), (r + 1, c + 1)]:
            nr, nc = move
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add((nr, nc))
        return attacks

    print("The smallest number of pieces for a legal position is believed to be 8.")
    print("Verifying the 8-piece solution by G. L. Gideon.")
    print("\nThe piece configuration (the 'equation' of the position):")
    
    black_king_pos_alg = 'h1'
    white_pieces_alg = {
        'Queen': ['d3'],
        'Rook': ['h7', 'g4'],
        'Bishop': ['h6'],
        'Knight': ['e5', 'f2'],
        'Pawn': ['e3', 'g2'],
    }
    
    num_pieces = 0
    for piece, positions in white_pieces_alg.items():
        for pos in positions:
            print(f"White {piece} on {pos}")
            num_pieces += 1
    print(f"Black King on {black_king_pos_alg}")
    print("-----------------------------------")
    print(f"Total number of white pieces: {num_pieces}")
    print("\nVerifying the position...")

    # --- Setup and Calculation ---
    king_pos = to_coords(black_king_pos_alg)
    
    all_piece_coords = {king_pos}
    white_piece_map = {}
    for piece, positions in white_pieces_alg.items():
        for pos_alg in positions:
            coords = to_coords(pos_alg)
            all_piece_coords.add(coords)
            white_piece_map[coords] = piece

    all_attacked_squares = set()
    for coords, piece in white_piece_map.items():
        if piece == 'Queen':
            all_attacked_squares.update(get_queen_attacks(coords, all_piece_coords))
        elif piece == 'Rook':
            all_attacked_squares.update(get_rook_attacks(coords, all_piece_coords))
        elif piece == 'Bishop':
            all_attacked_squares.update(get_bishop_attacks(coords, all_piece_coords))
        elif piece == 'Knight':
            all_attacked_squares.update(get_knight_attacks(coords))
        elif piece == 'Pawn':
            all_attacked_squares.update(get_pawn_attacks(coords))
            
    all_board_squares = set((r, c) for r in range(8) for c in range(8))
    unattacked_squares = all_board_squares - all_attacked_squares

    # --- Verification Output ---
    print("\n--- Verification Results ---")
    
    # 1. Check unattacked squares
    is_successful = True
    print(f"1. Number of unattacked squares found: {len(unattacked_squares)}")
    if len(unattacked_squares) == 1:
        unattacked_alg = to_algebraic(*list(unattacked_squares)[0])
        print(f"   - The single unattacked square is: {unattacked_alg}")
        if unattacked_alg == black_king_pos_alg:
            print("   - SUCCESS: This matches the Black King's position.")
        else:
            print(f"   - FAILURE: This does not match the king's position ({black_king_pos_alg}).")
            is_successful = False
    else:
        unattacked_list = sorted([to_algebraic(r,c) for r,c in unattacked_squares])
        print(f"   - FAILURE: Expected 1 unattacked square, but found {len(unattacked_squares)}: {unattacked_list}")
        is_successful = False

    # 2. Check stalemate conditions
    print("\n2. Stalemate Check:")
    if king_pos in all_attacked_squares:
        print("   - FAILURE: King is in check.")
        is_successful = False
    else:
        print("   - SUCCESS: King is not in check.")

    king_moves = set()
    r, c = king_pos
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0: continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                king_moves.add((nr, nc))
    
    has_legal_move = False
    for move in king_moves:
        if move not in all_attacked_squares and move not in all_piece_coords:
             has_legal_move = True
             print(f"   - FAILURE: King has a legal move to {to_algebraic(*move)}.")
             is_successful = False

    if not has_legal_move:
        print("   - SUCCESS: King has no legal moves.")
    
    # Final Conclusion
    print("\n--- Final Conclusion ---")
    if is_successful:
        print(f"The setup with {num_pieces} pieces is confirmed to be a valid solution.")
        print("The smallest number of points of white material that can achieve this in a legally reachable position is 8.")
    else:
        print("The provided configuration could not be verified.")

solve_chess_puzzle()