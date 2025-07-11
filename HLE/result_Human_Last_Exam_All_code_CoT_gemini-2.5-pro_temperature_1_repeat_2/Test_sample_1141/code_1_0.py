import numpy as np

def solve_cutting_problem():
    """
    Solves the billet cutting problem using a greedy heuristic.
    """
    # Billet dimensions in 0.5 cm grid units
    W, H, D = 32, 22, 8
    
    # Occupied space grid
    occupied = np.zeros((W, H, D), dtype=bool)

    # Item prices
    price_b2 = 150
    price_t1 = 5
    price_b1 = 1

    # Item counts and values
    n_b2, n_t1, n_b1 = 0, 0, 0
    value_b2, value_t1, value_b1 = 0, 0, 0

    # --- Step 1: Place B2 (radius 2cm = 4 units) ---
    b2_radius = 4
    b2_centers = [
        # Row 1
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        # Row 2
        (4, 12, 4), (12, 12, 4), (20, 12, 4), (28, 12, 4)
    ]
    
    n_b2 = len(b2_centers)
    value_b2 = n_b2 * price_b2
    
    # Mark occupied space for B2s
    for cx, cy, cz in b2_centers:
        # Define a bounding box for checking to improve performance
        min_x, max_x = max(0, cx - b2_radius), min(W, cx + b2_radius + 1)
        min_y, max_y = max(0, cy - b2_radius), min(H, cy + b2_radius + 1)
        min_z, max_z = max(0, cz - b2_radius), min(D, cz + b2_radius + 1)
        
        for i in range(min_x, max_x):
            for j in range(min_y, max_y):
                for k in range(min_z, max_z):
                    # Check if the center of the cell is within the sphere
                    if (i + 0.5 - cx)**2 + (j + 0.5 - cy)**2 + (k + 0.5 - cz)**2 <= b2_radius**2:
                        occupied[i, j, k] = True

    # --- Step 2: Place T1 (cube 1cm side = 2x2x2 units) ---
    # Center ranges for T1
    for cz in range(1, D - 1):
        for cy in range(1, H - 1):
            for cx in range(1, W - 1):
                # T1 occupies a 2x2x2 block of cells
                # The center (cx,cy,cz) requires cells [cx-1, cx] x [cy-1, cy] x [cz-1, cz]
                if not occupied[cx-1:cx+1, cy-1:cy+1, cz-1:cz+1].any():
                    occupied[cx-1:cx+1, cy-1:cy+1, cz-1:cz+1] = True
                    n_t1 += 1
    
    value_t1 = n_t1 * price_t1

    # --- Step 3: Place B1 (ball 1cm diameter = radius 0.5cm = 1 unit) ---
    # We approximate B1's occupied space as a 2x2x2 cube, same as T1
    # Center ranges for B1
    for cz in range(1, D - 1):
        for cy in range(1, H - 1):
            for cx in range(1, W - 1):
                 # B1 is approximated by a 2x2x2 block of cells
                 if not occupied[cx-1:cx+1, cy-1:cy+1, cz-1:cz+1].any():
                    occupied[cx-1:cx+1, cy-1:cy+1, cz-1:cz+1] = True
                    n_b1 += 1
    
    value_b1 = n_b1 * price_b1
    
    # --- Final Calculation ---
    total_value = value_b2 + value_t1 + value_b1
    
    print("Problem formulation is incorrect, but solving with a greedy algorithm.")
    print("Highest valid solution found:")
    print(f"B2 balls: {n_b2} x {price_b2} = {value_b2}")
    print(f"T1 cubes: {n_t1} x {price_t1} = {value_t1}")
    print(f"B1 balls: {n_b1} x {price_b1} = {value_b1}")
    print("\nFinal equation:")
    print(f"{n_b2}*150 + {n_t1}*5 + {n_b1}*1 = {total_value}")
    print(f"\nMaximum value found: <<< {total_value} >>>")

solve_cutting_problem()