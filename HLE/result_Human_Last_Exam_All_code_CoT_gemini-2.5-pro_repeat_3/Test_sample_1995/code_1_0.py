def get_attacked_squares(piece_type, row, col, occupied_squares):
    """
    Calculates the set of squares attacked by a piece.
    This simplified model assumes an empty board for attack paths.
    """
    attacked = set()
    
    # Rook moves (used by Rook and Queen)
    if piece_type in 'rq':
        for i in range(8):
            if i != row:
                attacked.add((i, col))
            if i != col:
                attacked.add((row, i))

    # Bishop moves (used by Bishop and Queen)
    if piece_type in 'bq':
        for i in range(1, 8):
            if 0 <= row + i < 8 and 0 <= col + i < 8: attacked.add((row + i, col + i))
            if 0 <= row + i < 8 and 0 <= col - i < 8: attacked.add((row + i, col - i))
            if 0 <= row - i < 8 and 0 <= col + i < 8: attacked.add((row - i, col + i))
            if 0 <= row - i < 8 and 0 <= col - i < 8: attacked.add((row - i, col - i))
            
    return attacked

def to_coords(alg_notation):
    """Converts algebraic notation (e.g., 'a1') to board coordinates (0, 0)."""
    file = ord(alg_notation[0]) - ord('a')
    rank = int(alg_notation[1]) - 1
    return rank, file

def to_alg(coords):
    """Converts board coordinates (0, 0) to algebraic notation 'a1'."""
    rank, file = coords
    return chr(ord('a') + file) + str(rank + 1)

def solve_chess_problem():
    """
    Verifies the solution to the chess problem and calculates the material cost.
    """
    # 1. Define the pieces, their positions, and their point values
    white_pieces = {
        'q': {'pos': 'd3', 'value': 9},
        'r': {'pos': 'e6', 'value': 5}
    }
    
    king_pos_alg = 'a1'
    king_pos_coords = to_coords(king_pos_alg)

    # Gather piece positions and calculate total material value
    occupied_coords = set()
    total_material = 0
    print("Proposed solution:")
    for piece_type, details in white_pieces.items():
        pos_alg = details['pos']
        occupied_coords.add(to_coords(pos_alg))
        total_material += details['value']
        print(f"- White {'Queen' if piece_type == 'q' else 'Rook'} on {pos_alg}")
    print(f"- Black King on safe square {king_pos_alg}\n")

    # 2. Verify stalemate condition
    print("Verifying stalemate condition for Black King on " + king_pos_alg + ":")
    king_adj_squares = {
        (r, c) for r in range(king_pos_coords[0] - 1, king_pos_coords[0] + 2)
        for c in range(king_pos_coords[1] - 1, king_pos_coords[1] + 2)
        if (0 <= r < 8 and 0 <= c < 8 and (r, c) != king_pos_coords)
    }
    
    all_attacked_by_white = set()
    for piece_type, details in white_pieces.items():
        r, c = to_coords(details['pos'])
        all_attacked_by_white.update(get_attacked_squares(piece_type, r, c, occupied_coords))
        
    stalemate = True
    for sq_coords in king_adj_squares:
        is_attacked = sq_coords in all_attacked_by_white
        print(f"- Square {to_alg(sq_coords)} adjacent to king is attacked: {is_attacked}")
        if not is_attacked:
            stalemate = False
    
    if stalemate:
        print("Result: King is stalemated.\n")
    else:
        print("Result: King is NOT stalemated.\n")

    # 3. Verify that all other 63 squares are controlled
    print("Verifying control of the board:")
    all_board_squares = {(r, c) for r in range(8) for c in range(8)}
    
    # Controlled squares are those attacked OR occupied by white pieces
    controlled_squares = all_attacked_by_white.union(occupied_coords)
    
    uncontrolled_squares = all_board_squares - controlled_squares
    
    print(f"Total squares on board: {len(all_board_squares)}")
    print(f"Squares controlled by White: {len(controlled_squares)}")
    print(f"Uncontrolled squares: {[to_alg(s) for s in uncontrolled_squares]}")

    if len(uncontrolled_squares) == 1 and king_pos_coords in uncontrolled_squares:
        print("Result: Verified. Exactly one square is uncontrolled, where the king is placed.\n")
    else:
        print("Result: Failed. The setup does not control exactly 63 squares.\n")

    # 4. Final calculation
    print("Calculating the smallest number of points:")
    equation_parts = [str(p['value']) for p in white_pieces.values()]
    print(f"{' + '.join(equation_parts)} = {total_material}")
    
    return total_material

if __name__ == '__main__':
    solve_chess_problem()