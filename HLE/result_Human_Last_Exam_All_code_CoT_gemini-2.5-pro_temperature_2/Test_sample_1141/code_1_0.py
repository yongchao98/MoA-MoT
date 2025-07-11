import sys

def solve_cutting_problem():
    """
    Solves the cutting problem based on the provided formulation.

    The analysis shows that B2 and T1 pieces are mutually exclusive.
    A B2-based solution (value > 1200) is far better than a T1-based solution (value ~ 20).
    Therefore, we calculate the number of B2 pieces and fill the rest of the space with B1 pieces.
    """

    # --- Piece Configuration ---
    # Prices
    val_b2 = 150
    val_b1 = 1
    val_t1 = 5
    
    # We can fit a maximum of 8 B2 balls in a 2x4 rectangular grid.
    b2_centers = [
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        (4, 12, 4), (12, 12, 4), (20, 12, 4), (28, 12, 4)
    ]
    num_b2 = len(b2_centers)
    num_t1 = 0 # Cannot be placed with B2 balls

    # --- Find Valid B1 Positions ---
    def is_b1_valid_vs_b2s(p, b2_list):
        """Check if a B1 position conflicts with any B2 ball."""
        x, y, z = p
        # Non-overlapping constraint for a B1 to a B2: (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= 25
        for xb, yb, zb in b2_list:
            if (x - xb)**2 + (y - yb)**2 + (z - zb)**2 < 25:
                return False
        return True

    # Generate all candidate points for B1 centers within the billet boundaries
    # that do not conflict with the pre-placed B2 balls.
    b1_candidates = []
    # B1/T1 center boundaries: 1 <= x <= 31, 1 <= y <= 21, 1 <= z <= 7
    for z in range(1, 8):
        for y in range(1, 22):
            for x in range(1, 32):
                p = (x, y, z)
                if is_b1_valid_vs_b2s(p, b2_centers):
                    b1_candidates.append(p)
    
    # --- Greedily Place B1 Balls ---
    b1_placed = []
    for p1 in b1_candidates:
        is_valid_placement = True
        x1, y1, z1 = p1
        for p2 in b1_placed:
            x2, y2, z2 = p2
            # Non-overlapping constraint for a B1 to a B1: (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= 4
            if (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2 < 4:
                is_valid_placement = False
                break
        if is_valid_placement:
            b1_placed.append(p1)

    num_b1 = len(b1_placed)
    
    # --- Calculate and Print Final Result ---
    total_value = (num_b2 * val_b2) + (num_b1 * val_b1) + (num_t1 * val_t1)
    
    print("Yes, the formulation can be solved as stated, though its constraints for T1 cubes are unusual.")
    print("The optimal strategy is to use B2 and B1 pieces only.")
    print(f"A greedy search found a solution with {num_b2} B2 balls and {num_b1} B1 balls.")
    print("The final calculation is:")
    print(f"{num_b2} * {val_b2} + {num_b1} * {val_b1} + {num_t1} * {val_t1} = {total_value}")
    
    # Return value for the final answer tag
    return total_value

if __name__ == '__main__':
    highest_value = solve_cutting_problem()
    # The required final answer tag.
    # To prevent it from being printed in the regular output, we use stderr.
    print(f"\n<<<ans>>>", file=sys.stderr)
    print(highest_value, file=sys.stderr)
    print(f"<<<ans>>>", file=sys.stderr)
