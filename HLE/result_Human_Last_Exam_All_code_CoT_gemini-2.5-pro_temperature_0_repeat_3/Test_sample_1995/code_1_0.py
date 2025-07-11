def get_queen_attacks(sq_r, sq_c):
    """Calculates all squares attacked by a queen."""
    attacks = set()
    # Rook moves
    for i in range(8):
        if i != sq_r: attacks.add((i, sq_c))
        if i != sq_c: attacks.add((sq_r, i))
    # Bishop moves
    for i in range(1, 8):
        if sq_r + i < 8 and sq_c + i < 8: attacks.add((sq_r + i, sq_c + i))
        else: break
    for i in range(1, 8):
        if sq_r + i < 8 and sq_c - i >= 0: attacks.add((sq_r + i, sq_c - i))
        else: break
    for i in range(1, 8):
        if sq_r - i >= 0 and sq_c + i < 8: attacks.add((sq_r - i, sq_c + i))
        else: break
    for i in range(1, 8):
        if sq_r - i >= 0 and sq_c - i >= 0: attacks.add((sq_r - i, sq_c - i))
        else: break
    return attacks

def get_pawn_attacks(sq_r, sq_c):
    """Calculates all squares attacked by a white pawn."""
    attacks = set()
    if sq_r < 7:
        if sq_c > 0: attacks.add((sq_r + 1, sq_c - 1))
        if sq_c < 7: attacks.add((sq_r + 1, sq_c + 1))
    return attacks

def square_to_an(r, c):
    """Converts (row, col) to algebraic notation like 'a1'."""
    return chr(ord('a') + c) + str(r + 1)

# --- Configuration ---
# White pieces
queen_pos = (5, 1)  # b6
pawn_pos = (6, 2)   # c7
white_pieces = {"Queen": queen_pos, "Pawn": pawn_pos}
white_occupied_squares = {queen_pos, pawn_pos}

# Black King
king_pos = (7, 0)   # a8

# --- Analysis ---
# 1. Calculate all attacks
queen_attacks = get_queen_attacks(queen_pos[0], queen_pos[1])
pawn_attacks = get_pawn_attacks(pawn_pos[0], pawn_pos[1])
all_white_attacks = queen_attacks | pawn_attacks

# 2. Check if the King is in check
king_is_in_check = king_pos in all_white_attacks
print(f"Position of Black King: {square_to_an(king_pos[0], king_pos[1])}")
print(f"Is the King in check? {king_is_in_check}")
print("-" * 20)

# 3. Check if all adjacent squares are attacked or occupied
kr, kc = king_pos
adjacent_squares = set()
for dr in [-1, 0, 1]:
    for dc in [-1, 0, 1]:
        if dr == 0 and dc == 0:
            continue
        nr, nc = kr + dr, kc + dc
        if 0 <= nr < 8 and 0 <= nc < 8:
            adjacent_squares.add((nr, nc))

print("Checking if King has any legal moves:")
has_legal_move = False
for r, c in adjacent_squares:
    is_attacked = (r, c) in all_white_attacks
    is_occupied = (r, c) in white_occupied_squares
    if not is_attacked and not is_occupied:
        has_legal_move = True
    print(f"  - Square {square_to_an(r, c)}: Attacked or Occupied = {is_attacked or is_occupied}")

print("-" * 20)
if not king_is_in_check and not has_legal_move:
    print("Result: The Black King is in stalemate.")
else:
    print("Result: The Black King is NOT in stalemate.")

# 4. Calculate the total point value
queen_points = 9
pawn_points = 1
total_points = queen_points + pawn_points
print("-" * 20)
print(f"The total point value of the white material is:")
print(f"{queen_points} (Queen) + {pawn_points} (Pawn) = {total_points}")
