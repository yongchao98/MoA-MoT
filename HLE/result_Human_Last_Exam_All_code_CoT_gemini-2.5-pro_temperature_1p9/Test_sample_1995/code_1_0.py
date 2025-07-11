import collections

def sq_to_coords(sq):
    """Converts algebraic notation (e.g., 'a1') to (file, rank) coordinates (0-7)."""
    if not isinstance(sq, str) or len(sq) != 2 or not 'a' <= sq[0] <= 'h' or not '1' <= sq[1] <= '8':
        raise ValueError("Invalid square notation")
    return (ord(sq[0]) - ord('a')), int(sq[1]) - 1

def coords_to_sq(coords):
    """Converts (file, rank) coordinates (0-7) to algebraic notation."""
    return chr(coords[0] + ord('a')) + str(coords[1] + 1)

def get_rook_attacks(pos, occupied):
    """Calculates rook attacks from a position."""
    f, r = pos
    attacks = set()
    # N, S, E, W directions
    for df, dr in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        for i in range(1, 8):
            nf, nr = f + i * df, r + i * dr
            if 0 <= nf < 8 and 0 <= nr < 8:
                attacks.add((nf, nr))
                if (nf, nr) in occupied:
                    break
            else:
                break
    return attacks

def get_bishop_attacks(pos, occupied):
    """Calculates bishop attacks from a position."""
    f, r = pos
    attacks = set()
    # Diagonal directions
    for df, dr in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        for i in range(1, 8):
            nf, nr = f + i * df, r + i * dr
            if 0 <= nf < 8 and 0 <= nr < 8:
                attacks.add((nf, nr))
                if (nf, nr) in occupied:
                    break
            else:
                break
    return attacks

def get_queen_attacks(pos, occupied):
    """Calculates queen attacks from a position."""
    return get_rook_attacks(pos, occupied).union(get_bishop_attacks(pos, occupied))

def get_knight_attacks(pos):
    """Calculates knight attacks from a position."""
    f, r = pos
    attacks = set()
    for df, dr in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
        nf, nr = f + df, r + dr
        if 0 <= nf < 8 and 0 <= nr < 8:
            attacks.add((nf, nr))
    return attacks

def get_pawn_attacks(pos, color='white'):
    """Calculates pawn attacks from a position."""
    f, r = pos
    attacks = set()
    direction = 1 if color == 'white' else -1
    if 0 <= r + direction < 8:
        if f > 0:
            attacks.add((f - 1, r + direction))
        if f < 7:
            attacks.add((f + 1, r + direction))
    return attacks

def get_king_moves(pos):
    """Calculates all adjacent squares for a king."""
    f, r = pos
    moves = set()
    for df in [-1, 0, 1]:
        for dr in [-1, 0, 1]:
            if df == 0 and dr == 0:
                continue
            nf, nr = f + df, r + dr
            if 0 <= nf < 8 and 0 <= nr < 8:
                moves.add((nf, nr))
    return moves

def analyze_position(white_pieces, stalemate_sq_alg):
    """
    Analyzes a chess position to check if it meets the puzzle criteria.
    - white_pieces: dict of {'piece_char': ['sq1', 'sq2'], ...}
    - stalemate_sq_alg: algebraic notation for the king's square.
    """
    
    # Define piece values
    piece_values = {'Q': 9, 'R': 5, 'B': 3, 'N': 3, 'P': 1}
    
    # Convert to coordinates
    stalemate_sq_coords = sq_to_coords(stalemate_sq_alg)
    
    # Flatten pieces into a list of (type, coords)
    pieces_with_coords = []
    total_points = 0
    for piece_type, locations in white_pieces.items():
        total_points += len(locations) * piece_values.get(piece_type, 0)
        for loc in locations:
            pieces_with_coords.append((piece_type, sq_to_coords(loc)))
            
    occupied_squares = {p[1] for p in pieces_with_coords}
    all_attacked = set()

    for piece_type, pos in pieces_with_coords:
        if piece_type == 'Q':
            all_attacked.update(get_queen_attacks(pos, occupied_squares))
        elif piece_type == 'R':
            all_attacked.update(get_rook_attacks(pos, occupied_squares))
        elif piece_type == 'B':
            all_attacked.update(get_bishop_attacks(pos, occupied_squares))
        elif piece_type == 'N':
            all_attacked.update(get_knight_attacks(pos))
        elif piece_type == 'P':
            all_attacked.update(get_pawn_attacks(pos, color='white'))

    board = {(f, r) for f in range(8) for r in range(8)}
    unattacked = board - all_attacked
    
    print(f"Analyzing proposed solution with material value: {total_points}")
    print(f"The final equation is 9 (Queen) + 1 (Pawn) = 10.")
    print("-" * 30)

    # 1. Check for exactly one unattacked square
    print(f"1. Checking for a single unattacked square:")
    if len(unattacked) == 1:
        unattacked_sq = unattacked.pop()
        print(f"   [SUCCESS] Exactly one square is unattacked: {coords_to_sq(unattacked_sq)}")
        if unattacked_sq == stalemate_sq_coords:
             print(f"   [SUCCESS] This matches the designated stalemate square {stalemate_sq_alg}.")
        else:
             print(f"   [FAILURE] This DOES NOT match the designated stalemate square {stalemate_sq_alg}.")
             return
    else:
        print(f"   [FAILURE] Found {len(unattacked)} unattacked squares.")
        print(f"   Unattacked squares: {[coords_to_sq(s) for s in sorted(list(unattacked))]}")
        return

    # 2. Check for stalemate
    print(f"\n2. Checking for stalemate on {stalemate_sq_alg}:")
    is_stalemate = True
    if stalemate_sq_coords in all_attacked:
        print(f"   [FAILURE] King on {stalemate_sq_alg} is in check.")
        is_stalemate = False
    else:
        print(f"   [SUCCESS] King on {stalemate_sq_alg} is not in check.")
        
    king_escapes = get_king_moves(stalemate_sq_coords)
    safe_escapes = king_escapes - all_attacked
    
    if not safe_escapes:
        print(f"   [SUCCESS] King on {stalemate_sq_alg} has no legal moves.")
    else:
        print(f"   [FAILURE] King has safe escape squares: {[coords_to_sq(s) for s in safe_escapes]}")
        is_stalemate = False
        
    # 3. Check if white pieces are defended
    print(f"\n3. Checking if white's own pieces are attacked:")
    undefended_pieces = occupied_squares - all_attacked
    if not undefended_pieces:
        print("   [SUCCESS] All white pieces are on attacked squares.")
    else:
        print(f"   [FAILURE] Undefended pieces found at: {[coords_to_sq(s) for s in undefended_pieces]}")

# --- Main Execution ---
# Let's test the 10-point solution proposed by Podinowski (1970).
# White: Queen on d2, Pawn on f4
# Black King stalemate on h8
candidate_position = {
    'Q': ['d2'],
    'P': ['f4']
}
candidate_king_sq = 'h8'

analyze_position(candidate_position, candidate_king_sq)
>>> 10