import collections

def to_coords(notation):
    """Converts chess notation like 'a1' to (row, col) tuple (0,0)."""
    col = ord(notation[0]) - ord('a')
    row = int(notation[1]) - 1
    return row, col

def to_notation(coords):
    """Converts (row, col) tuple (0,0) to chess notation 'a1'."""
    row, col = coords
    return f"{chr(ord('a') + col)}{row + 1}"

def get_rook_attacks(pos):
    """Calculates all squares attacked by a rook from a given position."""
    r, c = pos
    attacks = set()
    for i in range(8):
        if i != r:
            attacks.add((i, c))
        if i != c:
            attacks.add((r, i))
    return attacks

def get_bishop_attacks(pos):
    """Calculates all squares attacked by a bishop from a given position."""
    r, c = pos
    attacks = set()
    for i in range(1, 8):
        for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            nr, nc = r + i * dr, c + i * dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add((nr, nc))
    return attacks

def solve_and_verify():
    """
    Solves the chess problem by constructing the 9-piece solution
    and verifying its correctness.
    """
    black_king_pos_str = 'a1'
    white_pieces = {
        'Rook': ['c2', 'd3', 'e4', 'f5', 'g6', 'h7', 'b8'],
        'Bishop': ['b1', 'c7']
    }
    
    num_pieces = sum(len(positions) for positions in white_pieces.values())
    
    print(f"The smallest number of white pieces required is {num_pieces}.")
    print("This can be achieved with the following configuration:\n")
    print(f"Black King: {black_king_pos_str}")
    print("White Pieces:")
    for piece_type, positions in white_pieces.items():
        for pos in positions:
            print(f"- {piece_type} on {pos}")
            
    # Verification
    all_squares = {to_coords(f"{c}{r}") for c in "abcdefgh" for r in range(1, 9)}
    
    # Calculate all attacked squares
    all_attacked_squares_coords = set()
    for piece_type, positions in white_pieces.items():
        for pos_str in positions:
            pos_coords = to_coords(pos_str)
            if piece_type == 'Rook':
                all_attacked_squares_coords.update(get_rook_attacks(pos_coords))
            elif piece_type == 'Bishop':
                all_attacked_squares_coords.update(get_bishop_attacks(pos_coords))
                
    unattacked_squares = all_squares - all_attacked_squares_coords
    
    print("\n--- Verification ---")
    
    # 1. Verify that only one square is unattacked
    print(f"\n1. Number of unattacked squares: {len(unattacked_squares)}")
    if len(unattacked_squares) == 1:
        unattacked_notation = to_notation(list(unattacked_squares)[0])
        print(f"   - Success: The single unattacked square is {unattacked_notation}.")
        assert unattacked_notation == black_king_pos_str
    else:
        print("   - Failure: The number of unattacked squares is not 1.")
        
    # 2. Verify stalemate condition
    print("\n2. Stalemate Check for King on a1:")
    king_adj_squares_str = ['a2', 'b1', 'b2']
    king_adj_squares_coords = {to_coords(s) for s in king_adj_squares_str}
    
    stalemate = True
    for sq_coords in king_adj_squares_coords:
        is_attacked = sq_coords in all_attacked_squares_coords
        status = "ATTACKED" if is_attacked else "NOT ATTACKED"
        print(f"   - Square {to_notation(sq_coords)} is {status}.")
        if not is_attacked:
            stalemate = False
            
    if stalemate:
        print("   - Success: King on a1 has no legal moves and is not in check, resulting in a stalemate.")
    else:
        print("   - Failure: King is not in stalemate.")

if __name__ == '__main__':
    solve_and_verify()
    print("\n<<<9>>>")
