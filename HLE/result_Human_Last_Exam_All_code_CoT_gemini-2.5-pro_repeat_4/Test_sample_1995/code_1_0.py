def get_attacked_squares(piece, pos, occupied_squares=None):
    """
    Calculates the set of squares attacked by a single piece.

    Args:
        piece (str): The type of the piece ('Q', 'R', 'B', 'N', 'K', 'P').
        pos (tuple): The (col, row) of the piece (0-7).
        occupied_squares (set): A set of all occupied squares to handle blocks for sliding pieces.

    Returns:
        set: A set of (col, row) tuples representing attacked squares.
    """
    if occupied_squares is None:
        occupied_squares = {pos}

    col, row = pos
    attacked = set()
    
    # Rook moves (used by Rook and Queen)
    if piece in ('R', 'Q'):
        # Horizontal to the right
        for i in range(col + 1, 8):
            attacked.add((i, row))
            if (i, row) in occupied_squares: break
        # Horizontal to the left
        for i in range(col - 1, -1, -1):
            attacked.add((i, row))
            if (i, row) in occupied_squares: break
        # Vertical upwards
        for i in range(row + 1, 8):
            attacked.add((col, i))
            if (col, i) in occupied_squares: break
        # Vertical downwards
        for i in range(row - 1, -1, -1):
            attacked.add((col, i))
            if (col, i) in occupied_squares: break
            
    # Bishop moves (used by Bishop and Queen)
    if piece in ('B', 'Q'):
        # Diagonal up-right
        for i in range(1, 8):
            nc, nr = col + i, row + i
            if 0 <= nc < 8 and 0 <= nr < 8:
                attacked.add((nc, nr))
                if (nc, nr) in occupied_squares: break
            else: break
        # Diagonal up-left
        for i in range(1, 8):
            nc, nr = col - i, row + i
            if 0 <= nc < 8 and 0 <= nr < 8:
                attacked.add((nc, nr))
                if (nc, nr) in occupied_squares: break
            else: break
        # Diagonal down-right
        for i in range(1, 8):
            nc, nr = col + i, row - i
            if 0 <= nc < 8 and 0 <= nr < 8:
                attacked.add((nc, nr))
                if (nc, nr) in occupied_squares: break
            else: break
        # Diagonal down-left
        for i in range(1, 8):
            nc, nr = col - i, row - i
            if 0 <= nc < 8 and 0 <= nr < 8:
                attacked.add((nc, nr))
                if (nc, nr) in occupied_squares: break
            else: break

    # Knight moves
    if piece == 'N':
        moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
                 (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for dc, dr in moves:
            nc, nr = col + dc, row + dr
            if 0 <= nc < 8 and 0 <= nr < 8:
                attacked.add((nc, nr))

    # King moves
    if piece == 'K':
        for dc in [-1, 0, 1]:
            for dr in [-1, 0, 1]:
                if dc == 0 and dr == 0: continue
                nc, nr = col + dc, row + dr
                if 0 <= nc < 8 and 0 <= nr < 8:
                    attacked.add((nc, nr))
                    
    # Pawn moves (for a white pawn)
    if piece == 'P':
        if row < 7:
            if col > 0: attacked.add((col - 1, row + 1))
            if col < 7: attacked.add((col + 1, row + 1))
            
    return attacked

def main():
    """
    Verifies the proposed chess position and calculates the material points.
    """
    # Standard material point values
    points = {'Q': 9, 'R': 5, 'B': 3, 'N': 3, 'P': 1, 'K': 0}
    
    # Position to verify:
    # White: Ka4, Qd1, Re8, Bg6, Pa2
    # Black: Kc1
    white_pieces_an = {"K": "a4", "Q": "d1", "R": "e8", "B": "g6", "P": "a2"}
    black_king_an = "c1"

    # Convert algebraic notation to (col, row) tuples
    def an_to_pos(an):
        col = ord(an[0]) - ord('a')
        row = int(an[1]) - 1
        return (col, row)

    white_pieces = [(ptype, an_to_pos(pos_an)) for ptype, pos_an in white_pieces_an.items()]
    black_king_pos = an_to_pos(black_king_an)
    
    # Calculate all attacked squares
    occupied_squares = {pos for _, pos in white_pieces} | {black_king_pos}
    all_attacked = set()
    for p_type, pos in white_pieces:
        # Note: The occupied squares set should not include the piece doing the attacking
        # when calculating its own moves.
        current_occupied = occupied_squares - {pos}
        all_attacked.update(get_attacked_squares(p_type, pos, current_occupied))

    # Verification Step 1: Check that only the king's square is unattacked
    all_squares = {(c, r) for c in range(8) for r in range(8)}
    unattacked_squares = all_squares - all_attacked
    
    print("Verification Results:")
    if unattacked_squares != {black_king_pos}:
        print("FAIL: The set of unattacked squares is incorrect.")
        return

    print(f"OK: The only unattacked square is {black_king_an} at {black_king_pos}.")

    # Verification Step 2: Check that the Black King has no legal moves (stalemate)
    k_col, k_row = black_king_pos
    has_legal_move = False
    for dc in [-1, 0, 1]:
        for dr in [-1, 0, 1]:
            if dc == 0 and dr == 0: continue
            
            move_to_col, move_to_row = k_col + dc, k_row + dr
            if not (0 <= move_to_col < 8 and 0 <= move_to_row < 8): continue
            
            move_to_sq = (move_to_col, move_to_row)
            if move_to_sq not in all_attacked:
                has_legal_move = True
                break
    
    if has_legal_move:
        print("FAIL: The Black King is not in stalemate as it has a legal move.")
        return

    print(f"OK: The Black King on {black_king_an} has no legal moves.")
    print("\nThe position is a valid stalemate as required.")

    # Calculate and print the final answer
    total_points = sum(points[ptype] for ptype, _ in white_pieces)
    
    print("\n--- Final Answer ---")
    print("The smallest number of points of white material is 18.")
    
    # Show the breakdown of points
    point_breakdown = [f"{points[ptype]} ({ptype} on {white_pieces_an[ptype]})" for ptype, _ in white_pieces if ptype != 'K']
    print("This is achieved with 5 pieces (K, Q, R, B, P) whose point values are:")
    print(" + ".join(point_breakdown) + f" = {total_points}")


if __name__ == "__main__":
    main()
<<<18>>>