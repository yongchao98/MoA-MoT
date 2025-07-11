import math

def solve_billet_packing():
    """
    Solves the billet packing problem based on the user's formulation.

    The analysis shows that the provided non-overlapping constraints for T1 cubes
    are extremely restrictive. Specifically, the T1-B2 constraint makes it
    impossible for both piece types to exist in a solution. Therefore, the optimal
    strategy is to prioritize the high-value B2 pieces and fill the rest of the
    space with B1 pieces.

    This script implements a greedy algorithm to do that.
    """

    # --- Define Piece Properties and Billet Dimensions ---
    PIECES = {
        'B2': {'price': 150},
        'B1': {'price': 1},
        'T1': {'price': 5}
    }

    # --- Non-overlapping check function based on the user's exact formulation ---
    def check_no_overlap(p1, p2):
        """Checks if two pieces p1 and p2 can be placed. Returns True if they DON'T overlap."""
        pos1, t1 = p1['pos'], p1['type']
        pos2, t2 = p2['pos'], p2['type']

        dx = abs(pos1[0] - pos2[0])
        dy = abs(pos1[1] - pos2[1])
        dz = abs(pos1[2] - pos2[2])
        d2 = dx**2 + dy**2 + dz**2

        types = tuple(sorted((t1, t2)))

        if types == ('B1', 'B1'):
            return d2 >= 4
        if types == ('B1', 'B2'):
            return d2 >= 25
        if types == ('B2', 'B2'):
            return d2 >= 64
        # The T1 constraints are included for completeness but not used in the optimal strategy
        if types == ('B1', 'T1') or types == ('T1', 'T1'):
            return min(dx, dy, dz) >= 2
        if types == ('B2', 'T1'):
            return min(dx, dy, dz) >= 5
        return True

    # --- Main Greedy Algorithm ---
    placed_pieces = []
    total_value = 0
    
    # --- Phase 1: Place B2 pieces ---
    # A 4x2 grid of B2s can fit within the billet, separated by 8 units (2*radius)
    b2_centers = []
    for i in range(4):  # 4 B2s along the x-axis
        for j in range(2):  # 2 B2s along the y-axis
            x = 4 + i * 8
            y = 4 + j * 8
            z = 4  # B2s can only be at z=4
            b2_centers.append((x, y, z))

    num_b2 = 0
    for pos in b2_centers:
        piece = {'type': 'B2', 'price': PIECES['B2']['price'], 'pos': pos}
        placed_pieces.append(piece)
        total_value += piece['price']
        num_b2 += 1

    # --- Phase 2: Place B1 pieces in remaining space ---
    num_b1 = 0
    # Iterate through all possible B1 center locations
    for z in range(1, 8):
        for y in range(1, 22):
            for x in range(1, 32):
                new_b1 = {'type': 'B1', 'price': PIECES['B1']['price'], 'pos': (x, y, z)}
                is_valid = True
                
                # Check against all previously placed pieces
                for existing_piece in placed_pieces:
                    if not check_no_overlap(new_b1, existing_piece):
                        is_valid = False
                        break
                
                if is_valid:
                    placed_pieces.append(new_b1)
                    total_value += new_b1['price']
                    num_b1 += 1
    
    # --- Print Final Results ---
    print("Based on the provided formulation, the highest value is achieved by placing B2 and B1 pieces.")
    print(f"The final arrangement consists of {num_b2} B2 pieces and {num_b1} B1 pieces.")
    print("\nFinal Value Calculation:")
    # The final equation with each number explicitly shown
    print(f"{num_b2} * {PIECES['B2']['price']} (for B2) + {num_b1} * {PIECES['B1']['price']} (for B1) = {total_value}")


if __name__ == '__main__':
    solve_billet_packing()