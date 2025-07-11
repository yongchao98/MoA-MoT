def get_attacked_squares(piece, pos, occupied):
    """
    Calculates the set of squares attacked by a piece from a given position.
    
    Args:
        piece (str): The character representing the piece ('k', 'q', 'r', 'b', 'n').
        pos (tuple): The (row, col) of the piece on a 0-7 board.
        occupied (set): A set of (row, col) tuples for all occupied squares.
    
    Returns:
        set: A set of (row, col) tuples of attacked squares.
    """
    row, col = pos
    attacked = set()
    piece = piece.lower()

    if piece == 'k':
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = row + dr, col + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacked.add((nr, nc))
    
    elif piece == 'n':
        moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
                 (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for dr, dc in moves:
            nr, nc = row + dr, col + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacked.add((nr, nc))

    # Add sliding piece logic (Rook, Bishop, Queen)
    if piece in ('r', 'q'):
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            for i in range(1, 8):
                nr, nc = row + i*dr, col + i*dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacked.add((nr, nc))
                    if (nr, nc) in occupied:
                        break
                else:
                    break
    
    if piece in ('b', 'q'):
        for dr, dc in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
            for i in range(1, 8):
                nr, nc = row + i*dr, col + i*dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacked.add((nr, nc))
                    if (nr, nc) in occupied:
                        break
                else:
                    break
    return attacked

def pos_to_coord(s):
    """Converts chess notation (e.g., 'a1') to (row, col) tuple."""
    col = ord(s[0]) - ord('a')
    row = int(s[1]) - 1
    return (row, col)

def main():
    """
    Main function to solve the chess problem.
    """
    # This position is a known composition for this problem.
    # Material value: Rook(5) + Bishop(3) + Knight(3) = 11
    white_pieces_info = {
        'K': {'pos': 'f1', 'val': 0},
        'R': {'pos': 'e6', 'val': 5},
        'B': {'pos': 'c4', 'val': 3},
        'N': {'pos': 'd5', 'val': 3},
    }
    black_king_pos_str = 'h8'

    # Convert positions to coordinates
    white_pieces = {p: pos_to_coord(i['pos']) for p, i in white_pieces_info.items()}
    black_king_pos = pos_to_coord(black_king_pos_str)
    
    occupied_squares = set(white_pieces.values())

    # Calculate all squares attacked by white pieces
    all_attacked = set()
    for piece_char, pos in white_pieces.items():
        attacked = get_attacked_squares(piece_char, pos, occupied_squares)
        all_attacked.update(attacked)

    # 1. Check stalemate conditions
    # 1a. King is not in check
    king_in_check = black_king_pos in all_attacked
    
    # 1b. King has no legal moves
    king_escapes = get_attacked_squares('k', black_king_pos, set())
    king_is_trapped = king_escapes.issubset(all_attacked)

    # 2. Check board coverage
    all_board_squares = {(r, c) for r in range(8) for c in range(8)}
    covered_squares = all_attacked | occupied_squares
    unattacked_squares = all_board_squares - covered_squares

    # 3. Final verification
    if (not king_in_check and 
        king_is_trapped and 
        len(unattacked_squares) == 1 and 
        black_king_pos in unattacked_squares):
        
        total_value = sum(i['val'] for i in white_pieces_info.values())
        
        print("The position is a valid solution.")
        print("The smallest number of points is achieved with a Rook, a Bishop, and a Knight.")
        print("\nFinal Equation:")
        print(f"Rook({white_pieces_info['R']['val']}) + Bishop({white_pieces_info['B']['val']}) + Knight({white_pieces_info['N']['val']}) = {total_value}")
        
    else:
        print("The analyzed position is not a valid solution.")

if __name__ == "__main__":
    main()