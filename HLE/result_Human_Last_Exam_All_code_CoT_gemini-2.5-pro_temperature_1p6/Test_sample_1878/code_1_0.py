import time

def get_coords(sq):
    """Converts a square index (0-63) to (rank, file) tuple."""
    return sq // 8, sq % 8

def from_coords(r, f):
    """Converts (rank, file) tuple to a square index (0-63)."""
    if 0 <= r < 8 and 0 <= f < 8:
        return r * 8 + f
    return -1

def get_king_attacks(sq):
    """Returns a set of squares attacked by a king on the given square."""
    attacks = set()
    r, f = get_coords(sq)
    for dr in [-1, 0, 1]:
        for df in [-1, 0, 1]:
            if dr == 0 and df == 0: continue
            nr, nf = r + dr, f + df
            if 0 <= nr < 8 and 0 <= nf < 8:
                attacks.add(from_coords(nr, nf))
    return attacks

def get_knight_attacks(sq):
    """Returns a set of squares attacked by a knight on the given square."""
    attacks = set()
    r, f = get_coords(sq)
    for dr, df in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
        nr, nf = r + dr, f + df
        if 0 <= nr < 8 and 0 <= nf < 8:
            attacks.add(from_coords(nr, nf))
    return attacks

def get_pawn_attacks(sq):
    """Returns a set of squares attacked by a white pawn on the given square."""
    attacks = set()
    r, f = get_coords(sq)
    if r < 7:  # A pawn attacks forward to the next rank
        if f > 0: attacks.add(from_coords(r + 1, f - 1))
        if f < 7: attacks.add(from_coords(r + 1, f + 1))
    return attacks

def is_checkmate(k_sq, K_sq, N_sq, P_sq):
    """Checks if the position is checkmate against the black king."""
    # 1. Calculate all squares attacked by White
    K_attacks = get_king_attacks(K_sq)
    N_attacks = get_knight_attacks(N_sq)
    P_attacks = get_pawn_attacks(P_sq)
    white_attacks = K_attacks | N_attacks | P_attacks

    # 2. Confirm Black King is in check
    if k_sq not in N_attacks and k_sq not in P_attacks:
        return False

    # 3. Check if Black King has any legal moves
    for move_sq in get_king_attacks(k_sq):
        # A king can move to an empty square if it's not attacked
        if move_sq not in [K_sq, N_sq, P_sq] and move_sq not in white_attacks:
            return False
        # A king can capture a piece if the destination square is not re-attacked
        elif move_sq in [K_sq, N_sq, P_sq] and move_sq not in white_attacks:
            return False
            
    # If the king has no safe square to move to, it's checkmate.
    return True

def is_legal(k_sq, K_sq, N_sq, P_sq):
    """Checks if a checkmate position is legal by finding a valid prior move."""
    # Check if a move could have been made that was NOT from a position where Black was already in check.
    
    # 1. Retract Knight's move
    r_n, f_n = get_coords(N_sq)
    for dr, df in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
        prev_N_sq = from_coords(r_n - dr, f_n - df)
        if prev_N_sq != -1 and prev_N_sq not in [k_sq, K_sq, P_sq]:
            if not (k_sq in get_king_attacks(K_sq) or k_sq in get_pawn_attacks(P_sq)):
                if not (K_sq in get_king_attacks(k_sq)):
                    return True

    # 2. Retract King's move (must be a discovered check)
    for prev_K_sq in get_king_attacks(K_sq):
        if prev_K_sq not in [k_sq, N_sq, P_sq]:
            # Before the King moved, the Black King must not have been in check
            if not (k_sq in get_knight_attacks(N_sq) or k_sq in get_pawn_attacks(P_sq)):
                if not (prev_K_sq in get_king_attacks(k_sq)):
                    return True

    # 3. Retract Pawn's move (must be a non-capture)
    r_p, f_p = get_coords(P_sq)
    # One-step push
    prev_P_sq = from_coords(r_p - 1, f_p)
    if prev_P_sq != -1 and prev_P_sq not in [k_sq, K_sq, N_sq]:
        if not (k_sq in get_king_attacks(K_sq) or k_sq in get_knight_attacks(N_sq)):
            if not (K_sq in get_king_attacks(k_sq)):
                return True
    # Two-step push
    if r_p == 3: # Pawn is on rank 4
        prev_P_sq = from_coords(r_p - 2, f_p)
        intervening_sq = from_coords(r_p - 1, f_p)
        if prev_P_sq != -1 and prev_P_sq not in [k_sq, K_sq, N_sq] and intervening_sq not in [k_sq, K_sq, N_sq, P_sq]:
            if not (k_sq in get_king_attacks(K_sq) or k_sq in get_knight_attacks(N_sq)):
                if not (K_sq in get_king_attacks(k_sq)):
                    return True

    return False

def find_all_checkmates():
    """Finds all legal checkmate positions for KNP vs K."""
    legal_mates = set()
    # Optimization: A K+N mate requires the king to be on the edge. The pawn doesn't change this.
    edge_squares = {sq for sq in range(64) if get_coords(sq)[0] in [0, 7] or get_coords(sq)[1] in [0, 7]}

    print("Searching for legal checkmates... (This will take a few minutes)")
    start_time = time.time()
    
    for i, k_sq in enumerate(edge_squares):
        print(f"Checking Black King on edge square {i+1}/{len(edge_squares)}...", end='\r')
        for K_sq in range(64):
            if K_sq == k_sq or K_sq in get_king_attacks(k_sq): continue
            for N_sq in range(64):
                if N_sq == k_sq or N_sq == K_sq: continue
                # Pawn must be on ranks 2-7 (index 1-6)
                for r_p in range(1, 7):
                    for f_p in range(8):
                        P_sq = from_coords(r_p, f_p)
                        if P_sq in [k_sq, K_sq, N_sq]: continue
                        
                        if is_checkmate(k_sq, K_sq, N_sq, P_sq):
                            if is_legal(k_sq, K_sq, N_sq, P_sq):
                                norm_pos = tuple(sorted((k_sq, K_sq, N_sq, P_sq)))
                                legal_mates.add(norm_pos)
                                
    end_time = time.time()
    print(f"\nSearch complete. Time elapsed: {end_time - start_time:.2f} seconds.")
    return len(legal_mates)

if __name__ == '__main__':
    total_mates = find_all_checkmates()
    print(f"Total number of legal checkmates: {total_mates}")
