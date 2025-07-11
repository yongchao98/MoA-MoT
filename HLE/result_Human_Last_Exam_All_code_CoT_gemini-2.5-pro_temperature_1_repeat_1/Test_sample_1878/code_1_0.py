import itertools

def is_attacked_by_white(square_to_check, white_king_pos, white_knight_pos, white_pawn_pos):
    """
    Checks if a given square on the board is under attack by any of the white pieces.
    A square is also considered "attacked" by the king if the king is on an adjacent square.

    Args:
        square_to_check (tuple): The (row, col) of the square we're interested in.
        white_king_pos (tuple): The (row, col) of the White King.
        white_knight_pos (tuple): The (row, col) of the White Knight.
        white_pawn_pos (tuple): The (row, col) of the White Pawn.

    Returns:
        bool: True if the square is attacked, False otherwise.
    """
    r, c = square_to_check
    wk_r, wk_c = white_king_pos
    wn_r, wn_c = white_knight_pos
    wp_r, wp_c = white_pawn_pos

    # 1. Check for attack by the White King (adjacent squares)
    if max(abs(r - wk_r), abs(c - wk_c)) == 1:
        return True

    # 2. Check for attack by the White Knight ('L' shape move)
    dr_n = abs(r - wn_r)
    dc_n = abs(c - wn_c)
    if (dr_n == 1 and dc_n == 2) or (dr_n == 2 and dc_n == 1):
        return True

    # 3. Check for attack by the White Pawn
    # A white pawn at (wp_r, wp_c) attacks diagonally forward squares.
    # In our coordinate system (row 0 = 1st rank), this is the next row.
    if r == wp_r + 1 and abs(c - wp_c) == 1:
        return True

    return False

def count_checkmates():
    """
    Calculates the total number of legal checkmate positions with White (King, Knight, Pawn)
    against a lone Black King.
    """
    checkmate_count = 0
    squares = range(64)
    num_pieces = 4

    # Use itertools.permutations to efficiently generate all unique piece placements.
    for positions in itertools.permutations(squares, num_pieces):
        wk_sq, wn_sq, wp_sq, bk_sq = positions
        
        # Convert square indices (0-63) to (row, col) coordinates (0-7 for each)
        wk_pos = divmod(wk_sq, 8)
        wn_pos = divmod(wn_sq, 8)
        wp_pos = divmod(wp_sq, 8)
        bk_pos = divmod(bk_sq, 8)

        # === Legality Check 1: White pawn cannot be on the 1st or 8th rank. ===
        # (rows 0 and 7 in our 0-indexed coordinate system)
        if wp_pos[0] in (0, 7):
            continue

        # === Legality Check 2: Kings cannot be on adjacent squares. ===
        if max(abs(wk_pos[0] - bk_pos[0]), abs(wk_pos[1] - bk_pos[1])) <= 1:
            continue

        # === Checkmate Condition 1: The Black King must be in check. ===
        if not is_attacked_by_white(bk_pos, wk_pos, wn_pos, wp_pos):
            continue

        # === Checkmate Condition 2: The Black King must have no legal moves. ===
        has_safe_escape_square = False
        # Check all 8 potential squares around the Black King
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue  # This is the king's current square, not a move

                escape_r, escape_c = bk_pos[0] + dr, bk_pos[1] + dc

                # A move is legal if the destination is on the board AND not attacked by White.
                if 0 <= escape_r < 8 and 0 <= escape_c < 8:
                    escape_pos = (escape_r, escape_c)
                    if not is_attacked_by_white(escape_pos, wk_pos, wn_pos, wp_pos):
                        # Found a safe square to move to, so it's not checkmate.
                        has_safe_escape_square = True
                        break
            if has_safe_escape_square:
                break
        
        # If after checking all moves, no safe square was found, it is checkmate.
        if not has_safe_escape_square:
            checkmate_count += 1
            
    # The problem asks to output the numbers in the equation. Since the result is a single
    # number, we will print that number.
    print(checkmate_count)

if __name__ == '__main__':
    count_checkmates()