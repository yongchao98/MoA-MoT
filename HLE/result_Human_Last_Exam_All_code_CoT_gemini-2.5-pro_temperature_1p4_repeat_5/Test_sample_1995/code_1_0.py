def to_coords(sq):
    """Converts chess notation like 'a1' to (column, row) coordinates (0,0)."""
    if len(sq) != 2 or not ('a' <= sq[0] <= 'h') or not ('1' <= sq[1] <= '8'):
        raise ValueError("Invalid square notation")
    col = ord(sq[0]) - ord('a')
    row = int(sq[1]) - 1
    return col, row

def to_notation(coord):
    """Converts (column, row) coordinates (0,0) to chess notation 'a1'."""
    col, row = coord
    return f"{'abcdefgh'[col]}{row + 1}"

# This position was composed by G. N. Sachodjakin.
# Material cost: Queen (9) + Pawn (1) = 10 points.
white_pieces = {
    'K': to_coords('a2'),
    'Q': to_coords('b1'),
    'P': to_coords('e2')
}
black_king_pos = to_coords('a4')

print("Analyzing the position:")
print(f"White: King at a2, Queen at b1, Pawn at e2")
print(f"Black: King at a4 (the designated unattacked square)")
print("-" * 30)

# 1. Verify Stalemate
print("Verifying stalemate condition for Black King at a4:")
bk_escapes = {'a3', 'b3', 'b4', 'a5', 'b5'}
k_attacks_stalemate = {'a3', 'b3'}
q_attacks_stalemate = {'b4', 'a5', 'b5'}

# Check if all escape squares are covered.
if k_attacks_stalemate.union(q_attacks_stalemate) == bk_escapes:
    print("- The White King at a2 attacks a3 and b3.")
    print("- The White Queen at b1 attacks b4, a5, and b5.")
    print("- The square a4 itself is not attacked.")
    print("Result: The Black King is correctly stalemated.")
else:
    print("Error: The stalemate condition is not met.")
print("-" * 30)

# 2. Verify that all other 63 squares are attacked.
all_squares = {(c, r) for c in range(8) for r in range(8)}
occupied_by_white = set(white_pieces.values())
attacked_squares = set()

# King attacks
k_pos = white_pieces['K']
for dc in [-1, 0, 1]:
    for dr in [-1, 0, 1]:
        if dc == 0 and dr == 0: continue
        nc, nr = k_pos[0] + dc, k_pos[1] + dr
        if 0 <= nc < 8 and 0 <= nr < 8:
            attacked_squares.add((nc, nr))

# Pawn attacks
p_pos = white_pieces['P']
# A white pawn at (c,r) attacks (c-1, r+1) and (c+1, r+1)
for dc in [-1, 1]:
    nc, nr = p_pos[0] + dc, p_pos[1] + 1
    if 0 <= nc < 8 and 0 <= nr < 8:
        attacked_squares.add((nc, nr))

# Queen attacks (accounting for blocking by own pieces)
q_pos = white_pieces['Q']
directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1),
              (0, 1), (1, -1), (1, 0), (1, 1)]
for dc, dr in directions:
    for i in range(1, 8):
        nc, nr = q_pos[0] + dc * i, q_pos[1] + dr * i
        if not (0 <= nc < 8 and 0 <= nr < 8):
            break
        attacked_squares.add((nc, nr))
        # If the ray hits a white piece, it stops there.
        if (nc, nr) in occupied_by_white:
            break

# Determine unattacked squares
unattacked = sorted(list(all_squares - attacked_squares), key=lambda x: (x[1], x[0]))

print("Verifying that exactly 63 squares are attacked:")
print(f"Total squares attacked by White: {len(attacked_squares)}")
if len(unattacked) == 1:
    unattacked_sq_notation = to_notation(unattacked[0])
    print(f"The single unattacked square is: {unattacked_sq_notation}")
    if unattacked[0] == black_king_pos:
        print("Result: Verification successful. Exactly one square, a4, is unattacked.")
    else:
        print("Error: The unattacked square is not the Black King's square.")
else:
    print(f"Error: Found {len(unattacked)} unattacked squares instead of 1.")
    print(f"Unattacked squares: {[to_notation(sq) for sq in unattacked]}")

print("-" * 30)
print("Calculating the material cost:")
queen_val = 9
pawn_val = 1
total_val = queen_val + pawn_val
print(f"White Queen ({queen_val}) + White Pawn ({pawn_val}) = {total_val}")
print("\nThe smallest number of points of white material is 10.")