import math

def solve_packing():
    """
    Solves the packing problem for Qin Shi Huang's artifacts.
    """
    # Define problem parameters
    sheet_dims = [(140, 110), (110, 140)]  # Both orientations of the sheet
    square_dim = 10  # 10x10 cm square for the name
    circle_bbox_dim = 40  # 40x40 cm bounding box for the circle (20cm radius)

    square_chars = 4
    circle_chars = 999

    best_N = 0
    best_M = 0
    max_K = -1

    # Iterate through both sheet orientations
    for sheet_W, sheet_H in sheet_dims:
        # The best strategy is to maximize the number of circles (M) due to their high value.
        # We calculate the max number of circles that fit in a simple grid layout.
        
        # Max M for the current orientation
        max_m_cols = sheet_W // circle_bbox_dim
        max_m_rows = sheet_H // circle_bbox_dim
        M = max_m_cols * max_m_rows

        # Calculate the area occupied by this grid of circles
        occupied_W = max_m_cols * circle_bbox_dim
        occupied_H = max_m_rows * circle_bbox_dim

        # The remaining area is an L-shape if the circle grid is placed in a corner.
        # We can decompose this L-shape into two non-overlapping rectangles.
        # Rectangle 1
        rem_rect1_W = sheet_W
        rem_rect1_H = sheet_H - occupied_H
        # Rectangle 2
        rem_rect2_W = sheet_W - occupied_W
        rem_rect2_H = occupied_H

        # Calculate how many squares (N) fit in these remaining rectangles.
        n1 = (rem_rect1_W // square_dim) * (rem_rect1_H // square_dim)
        n2 = (rem_rect2_W // square_dim) * (rem_rect2_H // square_dim)
        N = n1 + n2

        # Calculate the total characters (K)
        K = (square_chars * N) + (circle_chars * M)

        # Update the best result found so far
        if K > max_K:
            max_K = K
            best_N = N
            best_M = M
            
    # Output the results
    print(f"Optimal number of squares (N): {best_N}")
    print(f"Optimal number of circles (M): {best_M}")
    print("\nThe maximum number of characters is calculated as follows:")
    # Output each number in the final equation
    print(f"{square_chars} * {best_N} + {circle_chars} * {best_M} = {max_K}")

    print("\nFinal Answer Format:")
    print(f"{best_N}:{best_M}:{max_K}")

# Run the solver
solve_packing()