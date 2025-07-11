import numpy as np

def solve_knot_tb():
    """
    Calculates the Thurston-Bennequin number for a knot defined by a grid diagram.
    """
    n = 5
    # O's are at (1,1), (2,2), (3,3), (4,4), (5,5)
    # X's are at (1,4), (2,5), (3,1), (4,2), (5,3)
    # Using 0-indexed coordinates (col, row)
    o_pos = [(i, i) for i in range(n)]
    x_pos = [(0, 3), (1, 4), (2, 0), (3, 1), (4, 2)]

    # --- 1. Calculate Writhe (w) ---
    o_col = {c: r for c, r in o_pos}
    x_col = {c: r for c, r in x_pos}
    o_row = {r: c for c, r in o_pos}
    x_row = {r: c for c, r in x_pos}

    writhe = 0
    for i in range(n):  # Column of vertical segment
        for j in range(n):  # Row of horizontal segment
            # Check for intersection
            y_min, y_max = min(o_col[i], x_col[i]), max(o_col[i], x_col[i])
            x_min, x_max = min(o_row[j], x_row[j]), max(o_row[j], x_row[j])
            
            if x_min < i < x_max and y_min < j < y_max:
                v_dir = np.sign(x_col[i] - o_col[i])
                h_dir = np.sign(x_row[j] - o_row[j])
                writhe += v_dir * h_dir

    # --- 2. Calculate Number of Components (c) ---
    # We trace the path: o --h--> x --v--> o
    # Map from an o's column to the next o's column
    path_map = {}
    for i in range(n):
        o_start_col, o_start_row = i, o_col[i]
        x_horiz_col, x_horiz_row = x_row[o_start_row], o_start_row
        o_end_col, o_end_row = x_horiz_col, o_col[x_horiz_col]
        path_map[o_start_col] = o_end_col

    num_components = 0
    visited = [False] * n
    for i in range(n):
        if not visited[i]:
            num_components += 1
            current = i
            while not visited[current]:
                visited[current] = True
                current = path_map[current]
    
    # --- 3. Calculate Rotation Number (rot) ---
    sum_y_o = sum(r for c, r in o_pos)
    sum_y_x = sum(r for c, r in x_pos)
    rot = (sum_y_x - sum_y_o) / 2

    # --- 4. Calculate Thurston-Bennequin Number (tb) ---
    tb = writhe - num_components - rot

    # --- Print the final equation ---
    print("Calculation Steps:")
    print(f"1. The writhe (w) is the sum of signs of all crossings.")
    print(f"   Calculated writhe w = {int(writhe)}")
    print(f"2. The number of components (c) is found by tracing the knot path.")
    print(f"   Calculated number of components c = {num_components}")
    print(f"3. The rotation number (rot) is calculated from the y-coordinates of the markings.")
    print(f"   Calculated rotation number rot = {int(rot)}")
    print("\nThe Thurston-Bennequin number (tb) is given by the formula: tb = w - c - rot")
    print("\nFinal Equation:")
    print(f"tb = {int(writhe)} - {num_components} - {int(rot)} = {int(tb)}")
    
    return int(tb)

# Run the solver and get the final answer
final_answer = solve_knot_tb()
print(f"\n<<<{final_answer}>>>")
