import sys

def to_coords(s):
    """Converts chess notation like 'a1' to a (row, col) tuple (0-indexed)."""
    if len(s) != 2 or not ('a' <= s[0] <= 'h') or not ('1' <= s[1] <= '8'):
        print(f"Invalid square notation: {s}")
        sys.exit(1)
    return (int(s[1:]) - 1, ord(s[0]) - ord('a'))

def to_square(r, c):
    """Converts a (row, col) tuple back to chess notation."""
    return chr(ord('a') + c) + str(r + 1)

def get_attack_set(piece_char, r_pos, c_pos):
    """
    Calculates the set of squares (row, col) attacked by a piece.
    This uses a simple "line of sight" definition, which is standard for such problems.
    """
    attacks = set()
    piece = piece_char.lower()
    # Rook and Queen straight-line moves
    if piece in ('r', 'q'):
        for i in range(8):
            if i != r_pos: attacks.add((i, c_pos))
            if i != c_pos: attacks.add((r_pos, i))
    # Bishop and Queen diagonal moves
    if piece in ('b', 'q'):
        for i in range(1, 8):
            if 0 <= r_pos + i < 8 and 0 <= c_pos + i < 8: attacks.add((r_pos + i, c_pos + i))
            if 0 <= r_pos + i < 8 and 0 <= c_pos - i < 8: attacks.add((r_pos + i, c_pos - i))
            if 0 <= r_pos - i < 8 and 0 <= c_pos + i < 8: attacks.add((r_pos - i, c_pos + i))
            if 0 <= r_pos - i < 8 and 0 <= c_pos - i < 8: attacks.add((r_pos - i, c_pos - i))
    # Knight moves
    elif piece == 'n':
        for dr, dc in [(1,2), (1,-2), (-1,2), (-1,-2), (2,1), (2,-1), (-2,1), (-2,-1)]:
            if 0 <= r_pos + dr < 8 and 0 <= c_pos + dc < 8: attacks.add((r_pos + dr, c_pos + dc))
    # King moves
    elif piece == 'k':
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                if 0 <= r_pos + dr < 8 and 0 <= c_pos + dc < 8: attacks.add((r_pos + dr, c_pos + dc))
    # Pawn (White) moves
    elif piece_char == 'P':
        if 0 <= r_pos + 1 < 8:
            if 0 <= c_pos - 1 < 8: attacks.add((r_pos + 1, c_pos - 1))
            if 0 <= c_pos + 1 < 8: attacks.add((r_pos + 1, c_pos + 1))
    return attacks

def solve_chess_puzzle():
    """
    Tests a candidate solution for the chess puzzle and prints the result.
    """
    # Candidate position: A self-protecting arrangement with a Queen and a Pawn.
    white_pieces = {'K': 'g7', 'Q': 'f6', 'P': 'c6'}
    black_king_square = 'a6'
    
    # Material cost calculation
    material_values = {'Q': 9, 'P': 1}
    material_cost = material_values['Q'] + material_values['P']

    # Convert positions to coordinate tuples for easier processing
    white_pieces_coords = {p: to_coords(s) for p, s in white_pieces.items()}
    black_king_coord = to_coords(black_king_square)

    # 1. Calculate the total set of attacked squares
    total_attacked_coords = set()
    for piece, pos_str in white_pieces.items():
        r, c = to_coords(pos_str)
        total_attacked_coords.update(get_attack_set(piece, r, c))

    # 2. Check for self-protection: All white pieces must be on attacked squares
    all_pieces_protected = True
    for piece_to_check, pos_coord in white_pieces_coords.items():
        # Get attacks from all *other* white pieces
        attackers_set = set()
        for other_piece, other_pos_str in white_pieces.items():
            if piece_to_check == other_piece: continue
            r, c = to_coords(other_pos_str)
            attackers_set.update(get_attack_set(other_piece, r, c))
        if pos_coord not in attackers_set:
            all_pieces_protected = False
            break

    # 3. Find the set of all unattacked squares
    all_board_coords = set((r, c) for r in range(8) for c in range(8))
    unattacked_coords = all_board_coords - total_attacked_coords

    # 4. Validate the final conditions
    # - Is there exactly one unattacked square?
    # - Is that square the Black King's square?
    # - Are all white pieces protected?
    is_solution_valid = (len(unattacked_coords) == 1 and
                         unattacked_coords.pop() == black_king_coord and
                         all_pieces_protected)

    # 5. Output the final result
    print("Verifying the solution...")
    if is_solution_valid:
        print("A valid position was found that meets all the problem's conditions.")
        print("The piece combination is a Queen and a Pawn.")
        # Output the numbers in the final equation as requested
        print(f"Queen ({material_values['Q']}) + Pawn ({material_values['P']}) = {material_cost}")
        print("\nThe smallest number of points of white material is:")
        print(material_cost)
    else:
        print("The proposed candidate solution did not work.")
        print("Further search would be required to find the correct minimal position.")

solve_chess_puzzle()