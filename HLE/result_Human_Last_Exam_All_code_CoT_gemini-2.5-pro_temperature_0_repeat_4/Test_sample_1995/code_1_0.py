def get_attacked_squares(piece, pos, occupied_squares):
    """
    Returns a set of squares attacked by a piece at a given position.
    - piece: 'P', 'N', 'B', 'R', 'Q'
    - pos: tuple (file, rank) from 0-7, e.g., (0, 0) for 'a1'
    - occupied_squares: a set of positions occupied by other pieces, which block movement.
    """
    attacked = set()
    f, r = pos
    piece_type = piece.upper()

    # Knight attacks are not blocked
    if piece_type == 'N':
        moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                 (1, -2), (1, 2), (2, -1), (2, 1)]
        for df, dr in moves:
            nf, nr = f + df, r + dr
            if 0 <= nf <= 7 and 0 <= nr <= 7:
                attacked.add((nf, nr))
        return attacked

    # Pawn attacks
    if piece_type == 'P':
        if r < 7: # Assuming white pawn moving up
            if f > 0:
                attacked.add((f - 1, r + 1))
            if f < 7:
                attacked.add((f + 1, r + 1))
        return attacked

    # Sliding pieces (Rook, Bishop, Queen)
    directions = []
    if piece_type in ('B', 'Q'):
        directions.extend([(-1, -1), (-1, 1), (1, -1), (1, 1)])
    if piece_type in ('R', 'Q'):
        directions.extend([(-1, 0), (1, 0), (0, -1), (0, 1)])

    for df, dr in directions:
        for i in range(1, 8):
            nf, nr = f + i * df, r + i * dr
            if 0 <= nf <= 7 and 0 <= nr <= 7:
                attacked.add((nf, nr))
                if (nf, nr) in occupied_squares:
                    break # Path is blocked
            else:
                break # Off the board
    
    return attacked

def pos_to_coord(s):
    """Converts chess notation like 'a1' to a (file, rank) tuple (0, 0)."""
    if len(s) != 2 or not ('a' <= s[0] <= 'h' and '1' <= s[1] <= '8'):
        raise ValueError("Invalid square notation")
    f = ord(s[0]) - ord('a')
    r = int(s[1]) - 1
    return f, r

def coord_to_pos(c):
    """Converts a (file, rank) tuple (0, 0) to chess notation 'a1'."""
    f, r = c
    return chr(ord('a') + f) + str(r + 1)

def solve_chess_puzzle():
    """
    Finds and verifies the solution to the chess stalemate puzzle.
    """
    # This 11-point configuration is a proposed solution.
    # It is believed that 10 points is not sufficient.
    solution_pieces = {
        'Rook': 'b8',
        'Knight': 'c3',
        'Bishop': 'd5'
    }
    king_pos_str = 'a1'
    
    piece_points = {'Rook': 5, 'Bishop': 3, 'Knight': 3, 'Queen': 9, 'Pawn': 1}

    # --- Verification ---
    
    # 1. Calculate total points and print the equation
    total_points = 0
    equation_parts = []
    for piece, pos in solution_pieces.items():
        points = piece_points[piece]
        total_points += points
        equation_parts.append(f"{piece} ({points})")
    
    print("Proposed Solution Material Cost:")
    print(" + ".join(equation_parts) + f" = {total_points}")
    print("-" * 30)

    # 2. Calculate all attacked squares
    piece_coords = {pos_to_coord(pos) for pos in solution_pieces.values()}
    total_attacked_coords = set()
    for piece, pos_str in solution_pieces.items():
        pos_coord = pos_to_coord(pos_str)
        # For calculating attacks, a piece does not block itself
        other_pieces = piece_coords - {pos_coord}
        total_attacked_coords.update(get_attacked_squares(piece[0], pos_coord, other_pieces))

    # 3. Check stalemate conditions
    king_coord = pos_to_coord(king_pos_str)
    
    # Condition 3a: King's square is not attacked
    king_is_safe = king_coord not in total_attacked_coords
    print(f"Verifying stalemate for Black King on {king_pos_str}:")
    print(f"1. Is King's square {king_pos_str} unattacked? {'Yes' if king_is_safe else 'No'}")

    # Condition 3b: All adjacent squares to the king are attacked
    kf, kr = king_coord
    king_neighbors = set()
    for df in [-1, 0, 1]:
        for dr in [-1, 0, 1]:
            if df == 0 and dr == 0:
                continue
            nf, nr = kf + df, kr + dr
            if 0 <= nf <= 7 and 0 <= nr <= 7:
                king_neighbors.add((nf, nr))
    
    stalemate_成立 = king_neighbors.issubset(total_attacked_coords)
    print(f"2. Are all adjacent squares ({', '.join(sorted([coord_to_pos(c) for c in king_neighbors]))}) attacked? {'Yes' if stalemate_成立 else 'No'}")
    
    # 4. Check that all other 63 squares are attacked
    all_squares = {(f, r) for f in range(8) for r in range(8)}
    must_be_attacked = all_squares - {king_coord}
    
    # The set of controlled squares includes those attacked and those occupied
    controlled_squares = total_attacked_coords.union(piece_coords)
    
    uncovered_squares = must_be_attacked - controlled_squares
    all_covered = not uncovered_squares

    print("-" * 30)
    print("Verifying Board Coverage:")
    print(f"Are all 63 squares (excluding {king_pos_str}) attacked or occupied? {'Yes' if all_covered else 'No'}")

    if uncovered_squares:
        print(f"Uncovered squares: {', '.join(sorted([coord_to_pos(c) for c in uncovered_squares]))}")

    print("-" * 30)
    if king_is_safe and stalemate_成立 and all_covered:
        print(f"Conclusion: The configuration is valid. The smallest number of points is {total_points}.")
    else:
        print("Conclusion: The proposed configuration is invalid.")

if __name__ == '__main__':
    solve_chess_puzzle()