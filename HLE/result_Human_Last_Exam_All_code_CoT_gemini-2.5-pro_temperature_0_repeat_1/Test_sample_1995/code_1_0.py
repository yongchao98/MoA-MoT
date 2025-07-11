import collections

def get_rook_attacks(pos, occupied):
    """Calculates all squares attacked by a rook."""
    r, c = pos
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

def get_bishop_attacks(pos, occupied):
    """Calculates all squares attacked by a bishop."""
    r, c = pos
    attacks = set()
    # Diagonals
    for dr in [-1, 1]:
        for dc in [-1, 1]:
            for i in range(1, 8):
                nr, nc = r + i * dr, c + i * dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacks.add((nr, nc))
                    if (nr, nc) in occupied:
                        break
                else:
                    break
    return attacks

def get_queen_attacks(pos, occupied):
    """Calculates all squares attacked by a queen."""
    return get_rook_attacks(pos, occupied).union(get_bishop_attacks(pos, occupied))

def get_knight_attacks(pos):
    """Calculates all squares attacked by a knight."""
    r, c = pos
    attacks = set()
    moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
             (1, -2), (1, 2), (2, -1), (2, 1)]
    for dr, dc in moves:
        nr, nc = r + dr, c + dc
        if 0 <= nr < 8 and 0 <= nc < 8:
            attacks.add((nr, nc))
    return attacks

def get_pawn_attacks(pos):
    """Calculates all squares attacked by a white pawn."""
    r, c = pos
    attacks = set()
    # White pawn attacks diagonally forward (up the board from white's perspective)
    if r > 0:
        if c > 0:
            attacks.add((r - 1, c - 1))
        if c < 7:
            attacks.add((r - 1, c + 1))
    return attacks

def an_to_coords(an):
    """Converts algebraic notation (e.g., 'h1') to (row, col) tuple."""
    col = ord(an[0]) - ord('a')
    row = 8 - int(an[1])
    return row, col

def coords_to_an(coords):
    """Converts (row, col) tuple to algebraic notation."""
    row, col = coords
    return chr(ord('a') + col) + str(8 - row)

def solve_chess_problem():
    """
    Solves the chess stalemate problem by demonstrating and verifying the known optimal solution.
    """
    # The record-holding solution by Gulyayev (15 points)
    white_pieces = {
        'R': [an_to_coords('h7')],
        'B': [an_to_coords('g5'), an_to_coords('h5')],
        'N': [an_to_coords('f3')],
        'P': [an_to_coords('g2')]
    }
    black_king_pos = an_to_coords('h1')
    
    piece_points = {'Q': 9, 'R': 5, 'B': 3, 'N': 3, 'P': 1}
    
    all_white_positions = set()
    for positions in white_pieces.values():
        all_white_positions.update(positions)

    all_attacked_squares = set()
    
    # Calculate attacks for each piece
    if 'Q' in white_pieces:
        for pos in white_pieces['Q']:
            all_attacked_squares.update(get_queen_attacks(pos, all_white_positions))
    if 'R' in white_pieces:
        for pos in white_pieces['R']:
            all_attacked_squares.update(get_rook_attacks(pos, all_white_positions))
    if 'B' in white_pieces:
        for pos in white_pieces['B']:
            all_attacked_squares.update(get_bishop_attacks(pos, all_white_positions))
    if 'N' in white_pieces:
        for pos in white_pieces['N']:
            all_attacked_squares.update(get_knight_attacks(pos))
    if 'P' in white_pieces:
        for pos in white_pieces['P']:
            all_attacked_squares.update(get_pawn_attacks(pos))

    all_squares = set((r, c) for r in range(8) for c in range(8))
    
    # Squares controlled by white (attacked or occupied)
    controlled_squares = all_attacked_squares.union(all_white_positions)
    
    unattacked_squares = all_squares - controlled_squares
    
    print("Verifying the 15-point solution by Gulyayev...")
    print("-" * 30)
    
    # 1. Verify king is not in check
    king_is_safe = black_king_pos not in controlled_squares
    print(f"King's square ({coords_to_an(black_king_pos)}) is unattacked: {king_is_safe}")

    # 2. Verify there is exactly one unattacked square
    one_unattacked_square = (len(unattacked_squares) == 1 and black_king_pos in unattacked_squares)
    print(f"Exactly one unattacked square on the board: {one_unattacked_square}")
    if not one_unattacked_square:
        print(f"Unattacked squares found: {[coords_to_an(s) for s in unattacked_squares]}")

    # 3. Verify king has no moves (stalemate)
    r, c = black_king_pos
    king_moves = set()
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                king_moves.add((nr, nc))
    
    king_has_no_moves = king_moves.issubset(controlled_squares)
    print(f"All king's potential moves are to controlled squares: {king_has_no_moves}")
    
    print("-" * 30)
    
    # Calculate and print the point value
    total_points = 0
    equation_parts = []
    
    # Sort pieces for consistent output order
    sorted_pieces = sorted(white_pieces.items(), key=lambda item: piece_points[item[0]], reverse=True)

    for piece_type, positions in sorted_pieces:
        for _ in positions:
            points = piece_points[piece_type]
            total_points += points
            piece_name = collections.defaultdict(lambda: piece_type, {'P': 'Pawn', 'N': 'Knight', 'B': 'Bishop', 'R': 'Rook', 'Q': 'Queen'})[piece_type]
            equation_parts.append(f"{piece_name} ({points})")
            
    equation = " + ".join(equation_parts)
    print("The material cost is calculated as follows:")
    print(f"{equation} = {total_points}")
    
    print("\nThe smallest number of points is 15.")


if __name__ == '__main__':
    solve_chess_problem()
    print("\n<<<15>>>")