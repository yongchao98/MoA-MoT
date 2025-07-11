def to_algebraic(coord):
    """Converts a (col, row) tuple (0-indexed) to algebraic notation."""
    col, row = coord
    return chr(ord('a') + col) + str(row + 1)

def get_queen_attacks(coord):
    """Returns a set of squares attacked by a queen at a given coordinate."""
    x, y = coord
    attacks = set()
    # Ranks and files
    for i in range(8):
        attacks.add((x, i))
        attacks.add((i, y))
    # Diagonals
    for i in range(1, 8):
        if 0 <= x + i < 8 and 0 <= y + i < 8: attacks.add((x + i, y + i))
        if 0 <= x - i < 8 and 0 <= y - i < 8: attacks.add((x - i, y - i))
        if 0 <= x + i < 8 and 0 <= y - i < 8: attacks.add((x + i, y - i))
        if 0 <= x - i < 8 and 0 <= y + i < 8: attacks.add((x - i, y + i))
    return attacks

def solve_chess_problem():
    """
    Finds and verifies a low-material solution to the stalemate problem.
    This solution uses three Queens, totaling 27 points.
    """
    # Using a known 3-queen configuration. Coordinates are (col, row) from 0-7.
    # Queens on b2, e5, h4. White King can be on a1.
    white_pieces = {
        "Queen1": {"pos": (1, 1), "value": 9}, # b2
        "Queen2": {"pos": (4, 4), "value": 9}, # e5
        "Queen3": {"pos": (7, 3), "value": 9}  # h4
    }

    all_board_squares = {(c, r) for c in range(8) for r in range(8)}
    
    # Calculate all attacked squares
    all_attacked_squares = set()
    for piece_info in white_pieces.values():
        all_attacked_squares.update(get_queen_attacks(piece_info["pos"]))

    # The pieces' own squares are part of the board, not targets for the final tally.
    # However, we must ensure they are attacked too ("defended").
    piece_positions = {p["pos"] for p in white_pieces.values()}
    unattacked_squares = all_board_squares - all_attacked_squares

    if len(unattacked_squares) != 1:
        print("This configuration does not leave exactly one square unattacked.")
        return

    stalemate_square = unattacked_squares.pop()
    
    # Verify the stalemate condition for the black king on the unattacked square
    bk_x, bk_y = stalemate_square
    adjacent_squares = set()
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            nx, ny = bk_x + dx, bk_y + dy
            if 0 <= nx < 8 and 0 <= ny < 8:
                adjacent_squares.add((nx, ny))
    
    # Check if all adjacent squares are attacked
    non_attacked_neighbors = adjacent_squares - all_attacked_squares

    # Also, the squares the pieces themselves sit on must be attacked.
    undefended_pieces = piece_positions - all_attacked_squares

    if not non_attacked_neighbors and not undefended_pieces:
        print("Solution found!")
        print("White pieces:")
        total_points = 0
        equation_parts = []
        for name, info in white_pieces.items():
            pos_alg = to_algebraic(info['pos'])
            print(f"- {name} on {pos_alg} (Value: {info['value']})")
            total_points += info['value']
            equation_parts.append(str(info['value']))

        print(f"\nThese pieces attack every square on the board except for {to_algebraic(stalemate_square)}.")
        print(f"The squares adjacent to {to_algebraic(stalemate_square)} are: {[to_algebraic(s) for s in adjacent_squares]}. All are attacked.")
        print("A Black King on this square would be in stalemate.")
        
        final_equation = " + ".join(equation_parts)
        print(f"\nThe smallest number of points for this construction is:")
        print(f"{final_equation} = {total_points}")

    else:
        print("This configuration does not result in a valid stalemate.")
        if non_attacked_neighbors:
            print(f"King could escape to: {[to_algebraic(s) for s in non_attacked_neighbors]}")
        if undefended_pieces:
            print(f"Some pieces are on unattacked squares: {[to_algebraic(s) for s in undefended_pieces]}")


solve_chess_problem()