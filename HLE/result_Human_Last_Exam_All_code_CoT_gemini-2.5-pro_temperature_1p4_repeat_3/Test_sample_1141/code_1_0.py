import math

def solve_cutting_problem():
    """
    Solves the cutting stock problem using a greedy algorithm based on the user's formulation.
    """
    
    # Piece definitions
    # Type 1: B1, radius 1, price 1
    # Type 2: B2, radius 4, price 150
    # Type 3: T1, half_side 1, price 5

    placed_pieces = []

    def check_overlap(p1, p2):
        """Checks if two pieces overlap based on the user's rules."""
        t1, c1 = p1['type'], p1['center']
        t2, c2 = p2['type'], p2['center']
        
        # Ensure t1 <= t2 for fewer cases in the conditional logic
        if t1 > t2:
            t1, t2 = t2, t1
            c1, c2 = c2, c1
            
        dx = abs(c1[0] - c2[0])
        dy = abs(c1[1] - c2[1])
        dz = abs(c1[2] - c2[2])
        
        d_sq = dx**2 + dy**2 + dz**2

        if t1 == 1 and t2 == 1: # B1 vs B1
            return d_sq < 4
        elif t1 == 1 and t2 == 2: # B1 vs B2
            return d_sq < 25
        elif t1 == 1 and t2 == 3: # B1 vs T1
            return min(dx, dy, dz) < 2
        elif t1 == 2 and t2 == 2: # B2 vs B2
            return d_sq < 64
        elif t1 == 2 and t2 == 3: # B2 vs T1
            return min(dx, dy, dz) < 5
        elif t1 == 3 and t2 == 3: # T1 vs T1
            return min(dx, dy, dz) < 2
        
        return False

    def can_place(new_piece, existing_pieces):
        """Checks if a new piece can be placed without overlapping existing ones."""
        for old_piece in existing_pieces:
            if check_overlap(new_piece, old_piece):
                return False
        return True

    # --- Phase 1: Place B2s (highest value) ---
    # B2 centers: 4 <= x <= 28, 4 <= y <= 18, z = 4
    for y in range(4, 18 + 1):
        for x in range(4, 28 + 1):
            z = 4
            new_piece = {'type': 2, 'center': (x, y, z), 'price': 150}
            if can_place(new_piece, placed_pieces):
                placed_pieces.append(new_piece)

    # --- Phase 2: Place T1s (next highest value-density) ---
    # T1 centers: 1 <= x <= 31, 1 <= y <= 21, 1 <= z <= 7
    for z in range(1, 7 + 1):
        for y in range(1, 21 + 1):
            for x in range(1, 31 + 1):
                new_piece = {'type': 3, 'center': (x, y, z), 'price': 5}
                if can_place(new_piece, placed_pieces):
                    placed_pieces.append(new_piece)
    
    # --- Phase 3: Place B1s ---
    # B1 centers: 1 <= x <= 31, 1 <= y <= 21, 1 <= z <= 7
    for z in range(1, 7 + 1):
        for y in range(1, 21 + 1):
            for x in range(1, 31 + 1):
                new_piece = {'type': 1, 'center': (x, y, z), 'price': 1}
                if can_place(new_piece, placed_pieces):
                    placed_pieces.append(new_piece)

    # --- Calculate and print results ---
    count_b2 = 0
    count_t1 = 0
    count_b1 = 0
    total_value = 0

    for piece in placed_pieces:
        total_value += piece['price']
        if piece['type'] == 2:
            count_b2 += 1
        elif piece['type'] == 3:
            count_t1 += 1
        elif piece['type'] == 1:
            count_b1 += 1
            
    print(f"{count_b2} * 150 + {count_b1} * 1 + {count_t1} * 5 = {total_value}")

solve_cutting_problem()