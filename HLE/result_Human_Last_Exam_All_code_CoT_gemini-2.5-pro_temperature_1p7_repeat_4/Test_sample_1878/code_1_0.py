import sys
import time

def solve():
    """
    This script calculates the number of unique, legal checkmate positions
    where White has a King, a Knight, and a Pawn, and Black has only a King.
    """
    
    # Step 1: Helper functions for board representation
    def to_coords(sq):
        """Converts a square index (0-63) to (row, col) coordinates."""
        return sq // 8, sq % 8

    def to_sq(r, c):
        """Converts (row, col) coordinates to a square index."""
        return r * 8 + c

    # Helper for symmetry transformations
    def transform_square(sq, t_idx):
        """Applies one of the 8 board symmetries to a square."""
        r, c = to_coords(sq)
        # Symmetries of a square: D4 group
        if t_idx == 0: tr, tc = r, c           # Identity
        elif t_idx == 1: tr, tc = c, 7 - r       # Rotate 90 degrees
        elif t_idx == 2: tr, tc = 7 - r, 7 - c   # Rotate 180 degrees
        elif t_idx == 3: tr, tc = 7 - c, r       # Rotate 270 degrees
        else:
            # Apply a vertical flip first, then rotations
            r, c = 7 - r, c
            if t_idx == 4: tr, tc = r, c           # Vertical Flip
            elif t_idx == 5: tr, tc = c, 7 - r       # Transpose
            elif t_idx == 6: tr, tc = 7 - r, 7 - c   # Horizontal Flip
            elif t_idx == 7: tr, tc = 7 - c, r       # Anti-Transpose
        return to_sq(tr, tc)

    # Pre-compute attack tables for all pieces for performance
    king_attacks = [set() for _ in range(64)]
    knight_attacks = [set() for _ in range(64)]
    pawn_attacks = [set() for _ in range(64)]  # White pawns move from high row to low row

    for sq in range(64):
        r, c = to_coords(sq)
        # King moves
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    king_attacks[sq].add(to_sq(nr, nc))
        # Knight moves
        for dr, dc in [(1,2), (1,-2), (-1,2), (-1,-2), (2,1), (2,-1), (-2,1), (-2,-1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                knight_attacks[sq].add(to_sq(nr, nc))
        # White Pawn attacks
        if r > 0:  # Pawn is not on the 8th rank (row 0)
            if c > 0: pawn_attacks[sq].add(to_sq(r - 1, c - 1))
            if c < 7: pawn_attacks[sq].add(to_sq(r - 1, c + 1))

    # Function to check if a square is attacked by White's army
    def is_attacked_by_white(target_sq, wk, wn, wp):
        if target_sq in king_attacks[wk]: return True
        if target_sq in knight_attacks[wn]: return True
        if target_sq in pawn_attacks[wp]: return True
        return False

    raw_mates = []
    
    # Step 2: Iterate through all possible piece placements
    for bk_pos in range(64):
        for wk_pos in range(64):
            if wk_pos == bk_pos: continue
            
            # Step 3: Legality Check 1 (King Proximity)
            if wk_pos in king_attacks[bk_pos]: continue

            for wn_pos in range(64):
                if wn_pos == bk_pos or wn_pos == wk_pos: continue

                # Pawns are only on ranks 2-7 (rows 1-6)
                for wp_row in range(1, 7):
                    for wp_col in range(8):
                        wp_pos = to_sq(wp_row, wp_col)
                        
                        if wp_pos in (bk_pos, wk_pos, wn_pos): continue

                        # Step 4: Check for Checkmate
                        # Condition 1: Black King must be in check
                        if not is_attacked_by_white(bk_pos, wk_pos, wn_pos, wp_pos):
                            continue
                        
                        # Condition 2: Black King must have no legal escape moves
                        can_escape = False
                        for move_sq in king_attacks[bk_pos]:
                            # Case A: Escape square is empty
                            if move_sq not in (wk_pos, wn_pos, wp_pos):
                                if not is_attacked_by_white(move_sq, wk_pos, wn_pos, wp_pos):
                                    can_escape = True
                                    break
                            # Case B: Escape square occupied by a capturable piece
                            else:
                                if move_sq == wn_pos: # Capture Knight?
                                    if not (wn_pos in king_attacks[wk_pos] or wn_pos in pawn_attacks[wp_pos]):
                                        can_escape = True # Can capture undefended knight
                                        break
                                elif move_sq == wp_pos: # Capture Pawn?
                                    if not (wp_pos in king_attacks[wk_pos] or wp_pos in knight_attacks[wn_pos]):
                                        can_escape = True # Can capture undefended pawn
                                        break
                        
                        if not can_escape:
                            raw_mates.append((bk_pos, wk_pos, wn_pos, wp_pos))

    # Step 5: Symmetry Reduction
    canonical_mates = set()
    for mate_pos in raw_mates:
        bk_pos, wk_pos, wn_pos, wp_pos = mate_pos
        
        symmetries = []
        for t in range(8):  # 8 transformations
            t_wp_pos = transform_square(wp_pos, t)
            
            # Transformed position is valid only if pawn remains on a legal rank
            if 1 <= to_coords(t_wp_pos)[0] <= 6:
                t_bk = transform_square(bk_pos, t)
                t_wk = transform_square(wk_pos, t)
                t_wn = transform_square(wn_pos, t)
                symmetries.append((t_bk, t_wk, t_wn, t_wp_pos))
        
        if symmetries:
            canonical_rep = min(symmetries)
            canonical_mates.add(canonical_rep)
    
    # Step 6: Final Output
    final_count = len(canonical_mates)
    num_raw_mates = len(raw_mates)
    
    # The prompt asks for an equation, which can be interpreted as showing the work.
    print(f"Number of raw checkmate positions found: {num_raw_mates}")
    print(f"After symmetry reduction, the final equation is:")
    print(f"Total number of unique checkmates = {final_count}")
    
    return final_count

# Execute the solver
final_answer = solve()
print(f"<<<{final_answer}>>>")