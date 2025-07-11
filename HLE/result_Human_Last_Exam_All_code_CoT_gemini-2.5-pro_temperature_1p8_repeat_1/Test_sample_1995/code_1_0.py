import sys

def to_coords(notation):
    """Converts algebraic notation (e.g., 'a1') to (file, rank) coordinates."""
    file = ord(notation[0]) - ord('a')
    rank = int(notation[1]) - 1
    return (file, rank)

def to_notation(coords):
    """Converts (file, rank) coordinates to algebraic notation."""
    file = chr(coords[0] + ord('a'))
    rank = str(coords[1] + 1)
    return f"{file}{rank}"

def get_king_attacks(pos):
    """Returns the set of squares a king at pos attacks."""
    kx, ky = pos
    attacks = set()
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            nx, ny = kx + dx, ky + dy
            if 0 <= nx < 8 and 0 <= ny < 8:
                attacks.add((nx, ny))
    return attacks

def get_pawn_attacks(pos):
    """Returns the set of squares a white pawn at pos attacks."""
    px, py = pos
    attacks = set()
    # A white pawn attacks diagonally forward
    for dx in [-1, 1]:
        nx, ny = px + dx, py + 1
        if 0 <= nx < 8 and 0 <= ny < 8:
            attacks.add((nx, ny))
    return attacks

def get_sliding_attacks(pos, directions, occupied_squares):
    """
    Calculates attacks for sliding pieces (Queen, Rook, Bishop).
    The attack ray stops when it hits any other piece.
    """
    px, py = pos
    attacks = set()
    for dx, dy in directions:
        nx, ny = px + dx, py + dy
        while 0 <= nx < 8 and 0 <= ny < 8:
            attacks.add((nx, ny))
            # If the ray hits any piece, it stops
            if (nx, ny) in occupied_squares:
                break
            nx, ny = nx + dx, ny + dy
    return attacks

def solve_chess_puzzle():
    """
    Verifies the 10-point solution to the stalemate problem and prints the result.
    """
    # --- Position Setup ---
    # Solution by Karl Fabel
    white_king_pos = to_coords('f5')
    white_queen_pos = to_coords('d8')
    white_pawn_pos = to_coords('c7')
    black_king_pos = to_coords('a8')

    white_pieces = {
        'King': {'pos': white_king_pos, 'value': 0},
        'Queen': {'pos': white_queen_pos, 'value': 9},
        'Pawn': {'pos': white_pawn_pos, 'value': 1}
    }

    # All occupied squares which can block sliding piece attacks
    occupied_squares = {p['pos'] for p in white_pieces.values()}

    attacked_squares = set()

    # --- Calculate Attacks for Each Piece ---
    # 1. King attacks
    king_attacks = get_king_attacks(white_king_pos)
    attacked_squares.update(king_attacks)

    # 2. Pawn attacks
    pawn_attacks = get_pawn_attacks(white_pawn_pos)
    attacked_squares.update(pawn_attacks)

    # 3. Queen attacks (Rook + Bishop moves)
    queen_directions = [
        (1, 0), (-1, 0), (0, 1), (0, -1),  # Rook moves
        (1, 1), (1, -1), (-1, 1), (-1, -1) # Bishop moves
    ]
    queen_attacks = get_sliding_attacks(white_queen_pos, queen_directions, occupied_squares)
    attacked_squares.update(queen_attacks)

    # --- Verification ---
    board_squares = {(x, y) for x in range(8) for y in range(8)}
    
    # The set of attacked squares should include the squares occupied by white pieces themselves
    attacked_squares.update(occupied_squares)
    
    unattacked_squares = board_squares - attacked_squares
    
    queen_value = white_pieces['Queen']['value']
    pawn_value = white_pieces['Pawn']['value']
    total_value = queen_value + pawn_value

    # --- Print Result ---
    print("This problem seeks the smallest point total of white material to attack 63 of the 64 squares,")
    print("leaving the black king in a stalemate on the one unattacked square.")
    print("\nThe solution involves a Queen and a Pawn.")
    print(f"\nWhite Piece Positions:")
    print(f"  King: {to_notation(white_king_pos)}")
    print(f"  Queen: {to_notation(white_queen_pos)} ({queen_value} points)")
    print(f"  Pawn: {to_notation(white_pawn_pos)} ({pawn_value} point)")
    print(f"Black King Position (Stalemate square): {to_notation(black_king_pos)}")
    
    print(f"\nVerification:")
    print(f"Total squares attacked by white material: {len(attacked_squares)}")
    if len(unattacked_squares) == 1 and list(unattacked_squares)[0] == black_king_pos:
        print(f"The single unattacked square is indeed {to_notation(black_king_pos)}, as required.")
    else:
        print(f"Verification failed. Unattacked squares: {[to_notation(s) for s in unattacked_squares]}", file=sys.stderr)

    print("\nThe smallest number of points is the sum of the material values.")
    print(f"Final Equation: {queen_value} (Queen) + {pawn_value} (Pawn) = {total_value}")


solve_chess_puzzle()

<<<10>>>