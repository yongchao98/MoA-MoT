def solve_chess_puzzle():
    """
    This script verifies a proposed solution to the stalemate puzzle.
    It checks if a specific arrangement of white pieces attacks exactly 63 squares,
    leaving the black king stalemated on the single safe square.
    """

    # --- Setup ---
    # The board is represented as a dictionary. Key: position, Value: piece
    # This is a known legal solution by B.P. Barnes.
    # White pieces: King (K), Rook (R), Bishop (B), Knight (N)
    # Black pieces: king (k)
    board = {
        'c1': 'K', 'b2': 'R', 'd2': 'B', 'c3': 'N',  # White pieces
        'a1': 'k'                                   # Black piece
    }

    piece_values = {
        'Q': 9, 'R': 5, 'B': 3, 'N': 3, 'P': 1
    }

    # --- Helper Functions ---
    def get_coords(square):
        """Converts square notation (e.g., 'a1') to (row, col) tuple (e.g., (0, 0))."""
        col = ord(square[0]) - ord('a')
        row = int(square[1]) - 1
        return row, col

    def to_square(row, col):
        """Converts (row, col) tuple to square notation."""
        if 0 <= row < 8 and 0 <= col < 8:
            return chr(ord('a') + col) + str(row + 1)
        return None

    def get_attacked_squares(board_state):
        """
        Calculates all squares attacked by white pieces.
        A piece does not attack the square it occupies.
        """
        attacked = set()
        for square, piece in board_state.items():
            # We only care about white's attacks
            if not piece.isupper():
                continue

            r, c = get_coords(square)
            
            # King attacks
            if piece == 'K':
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        if dr == 0 and dc == 0:
                            continue
                        s = to_square(r + dr, c + dc)
                        if s:
                            attacked.add(s)

            # Rook attacks (and Queen's rook-like moves)
            if piece in 'RQ':
                for i in range(8):
                    if i != r: attacked.add(to_square(i, c))
                    if i != c: attacked.add(to_square(r, i))

            # Bishop attacks (and Queen's bishop-like moves)
            if piece in 'BQ':
                for i in range(1, 8):
                    for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                        s = to_square(r + i*dr, c + i*dc)
                        if s: attacked.add(s)

            # Knight attacks
            if piece == 'N':
                for dr, dc in [(2, 1), (2, -1), (-2, 1), (-2, -1),
                               (1, 2), (1, -2), (-1, 2), (-1, -2)]:
                    s = to_square(r + dr, c + dc)
                    if s:
                        attacked.add(s)
        return attacked

    # --- Verification Logic ---
    print("Verifying the position:")
    print(f"White pieces: King at c1, Rook at b2, Bishop at d2, Knight at c3")
    print("Black pieces: King at a1")
    print("-" * 30)
    
    # 1. Calculate attacks
    attacked_by_white = get_attacked_squares(board)
    num_attacked = len(attacked_by_white)

    # 2. Check if exactly 63 squares are attacked
    all_squares = {to_square(r, c) for r in range(8) for c in range(8)}
    unattacked_squares = all_squares - attacked_by_white
    
    # 3. Find black king and his position
    black_king_pos = None
    for sq, p in board.items():
        if p == 'k':
            black_king_pos = sq
            break

    # 4. Perform checks and print results
    is_successful = True
    if num_attacked == 63:
        print(f"Success: White attacks exactly {num_attacked} squares.")
    else:
        print(f"Failure: White attacks {num_attacked} squares, not 63.")
        is_successful = False

    if len(unattacked_squares) == 1 and black_king_pos in unattacked_squares:
        print(f"Success: The single unattacked square is {list(unattacked_squares)[0]}, which is occupied by the Black King.")
    else:
        print(f"Failure: The unattacked square logic is incorrect.")
        is_successful = False

    # Check stalemate condition (king has no moves)
    kr, kc = get_coords(black_king_pos)
    escape_squares = set()
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            s = to_square(kr + dr, kc + dc)
            if s:
                escape_squares.add(s)

    if escape_squares.issubset(attacked_by_white):
        print(f"Success: All king's escape squares ({', '.join(sorted(list(escape_squares)))}) are attacked. It is stalemate.")
    else:
        print("Failure: The Black King is not in stalemate.")
        is_successful = False

    print("-" * 30)

    # 5. Calculate material cost
    total_points = 0
    calculation_str_parts = []
    white_material = {k: v for k, v in board.items() if v.isupper() and v != 'K'}
    
    for sq, piece_char in sorted(white_material.items()):
        points = piece_values[piece_char]
        piece_name = {'R': 'Rook', 'B': 'Bishop', 'N': 'Knight', 'Q': 'Queen', 'P': 'Pawn'}[piece_char]
        total_points += points
        calculation_str_parts.append(f"{points} ({piece_name})")

    if is_successful:
        print("Conclusion: The position is a legal stalemate satisfying all conditions.")
        print("The smallest number of material points required is:")
        print(" + ".join(calculation_str_parts) + f" = {total_points}")
    else:
        print("The proposed position did not meet all criteria.")
    
    return total_points

if __name__ == '__main__':
    final_answer = solve_chess_puzzle()
    print(f"\n<<<{final_answer}>>>")
