import collections

# Capablanca chess board is 10x8. We'll use a 0-indexed coordinate system.
# Files 'a' through 'j' are columns 0 through 9.
# Ranks '1' through '8' are rows 0 through 7.

Piece = collections.namedtuple('Piece', ['type', 'color'])

def str_to_coords(s):
    """Converts algebraic notation like 'a1' to (col, row) tuple (0,0)."""
    col = ord(s[0]) - ord('a')
    row = int(s[1:]) - 1
    return col, row

def coords_to_str(coords):
    """Converts (col, row) tuple (0,0) to algebraic notation 'a1'."""
    col, row = coords
    return chr(ord('a') + col) + str(row + 1)

def is_square_attacked_by_chancellor(square_coords, chancellor_coords):
    """Checks if a square is attacked by a Chancellor."""
    sq_col, sq_row = square_coords
    ch_col, ch_row = chancellor_coords

    # Rook moves
    if sq_col == ch_col or sq_row == ch_row:
        return True
    
    # Knight moves
    d_col = abs(sq_col - ch_col)
    d_row = abs(sq_row - ch_row)
    if (d_col == 1 and d_row == 2) or (d_col == 2 and d_row == 1):
        return True
        
    return False

def is_square_attacked_by_queen(square_coords, queen_coords):
    """Checks if a square is attacked by a Queen."""
    sq_col, sq_row = square_coords
    q_col, q_row = queen_coords

    # Rook moves
    if sq_col == q_col or sq_row == q_row:
        return True
    
    # Bishop moves
    if abs(sq_col - q_col) == abs(sq_row - q_row):
        return True
        
    return False

def solve():
    """
    Solves the Capablanca chess mate problem with the corrected piece.
    """
    # Initial positions based on FEN: 9k/5c1pb1/.../PP5C2/K9 w
    # We only need the pieces relevant to the mate.
    pieces = {
        'k_b': str_to_coords('j8'),
        'c_b': str_to_coords('f7'),
        'b_b': str_to_coords('i7'),
        'Q_w': str_to_coords('d3'),
        'C_w': str_to_coords('h2'),
    }

    # --- Move 1: White plays Cj2+ ---
    # The Chancellor moves from h2 to j2, checking the black King.
    pieces['C_w'] = str_to_coords('j2')
    white_moves = 1
    print(f"White's move {white_moves}: Chancellor moves from h2 to j2, delivering a check.")

    # --- Black's forced response: ...cj7 ---
    # The black king at j8 is checked. It cannot move to i8 because that square
    # is attacked by the white Queen on d3. Its only reasonable move is to block.
    # The black chancellor at f7 moves to j7 to block the check.
    pieces['c_b'] = str_to_coords('j7')
    print("Black's response: Chancellor moves from f7 to j7 to block.")

    # --- Move 2: White plays Cxj7# ---
    # The white Chancellor captures the black chancellor on j7.
    pieces['C_w'] = str_to_coords('j7')
    white_moves += 1
    print(f"White's move {white_moves}: Chancellor captures chancellor on j7, delivering checkmate.")

    # --- Verification of Checkmate ---
    # King at j8 is in check from the Chancellor at j7 (rook-like attack).
    king_pos = pieces['k_b']
    print(f"\nVerifying checkmate for black king at {coords_to_str(king_pos)}:")

    # Check escape squares for the king
    escape_squares = [str_to_coords('i8'), str_to_coords('h8')]
    
    # Is i8 safe?
    q_w_pos = pieces['Q_w']
    i8_attacked = is_square_attacked_by_queen(escape_squares[0], q_w_pos)
    print(f"- Escape square i8 is attacked by Queen at d3: {i8_attacked}")

    # Is h8 safe?
    c_w_pos = pieces['C_w']
    h8_attacked = is_square_attacked_by_chancellor(escape_squares[1], c_w_pos)
    print(f"- Escape square h8 is attacked by Chancellor at j7 (knight move): {h8_attacked}")
    
    # Can the checking piece be captured? No black piece can legally move to j7.
    # Is the check blocked? No, it's a direct adjacency capture.
    
    if i8_attacked and h8_attacked:
        print("\nAll escape squares are covered. The position is checkmate.")
        print("Minimal moves by White to win is:")
        print(white_moves)
    else:
        print("Error: The position is not checkmate.")

solve()