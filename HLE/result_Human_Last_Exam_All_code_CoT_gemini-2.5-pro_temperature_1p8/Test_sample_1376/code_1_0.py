import math

def calculate_max_artifacts():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # Define constants based on the problem description
    sheet_w = 140
    sheet_h = 110
    square_side = 10
    circle_radius = 20
    circle_diameter = 40
    chars_per_square = 4
    chars_per_circle = 999

    print("--- Step 1 & 2: Define Goal and Strategy ---")
    print(f"The goal is to maximize the total characters K = {chars_per_square}*N + {chars_per_circle}*M.")
    print("Strategy: Prioritize circles (M) due to their high character count.")
    print("\n--- Step 3: Calculate the Maximum Number of Circles (M) ---")

    def get_max_circles(w, h, r):
        """Calculates max circles with hexagonal packing."""
        d = 2 * r
        # Circles in a non-staggered row
        n1 = math.floor(w / d)
        # Circles in a staggered row
        n2 = math.floor((w - r) / d)
        # Number of rows
        if h < d:
            return 0
        m = 1 + math.floor((h - d) / (r * math.sqrt(3)))
        
        # Total circles depending on row parity
        if m == 0:
            return 0
        elif m % 2 == 1:
            return (m + 1) // 2 * n1 + (m - 1) // 2 * n2
        else: # m is even
            # We can start with the shorter or longer row type
            c1 = (m / 2) * n1 + (m / 2) * n2 # starting with n1
            c2 = (m / 2) * n2 + (m / 2) * n1 # starting with n2 (same result)
            return c1


    # Orientation 1: 140x110
    m_orient1 = get_max_circles(sheet_w, sheet_h, circle_radius)
    print(f"With sheet orientation {sheet_w}x{sheet_h}, we can fit {m_orient1} circles.")

    # Orientation 2: 110x140
    m_orient2 = get_max_circles(sheet_h, sheet_w, circle_radius)
    print(f"With sheet orientation {sheet_h}x{sheet_w}, we can fit {m_orient2} circles.")
    
    M = max(m_orient1, m_orient2)
    print(f"The maximum number of circles (M) is {M}.")

    print("\n--- Step 4: Calculate the Number of Squares (N) in Leftover Space ---")
    # This calculation assumes the optimal M=9 configuration (3 staggered rows of 3)
    # The sheet is oriented as 140x110.
    N = 0
    if M == 9:
        print("For M=9, circles are arranged in 3 staggered rows along the 140cm side.")
        # Calculate vertical positions of rows
        y_dist_rows = circle_radius * math.sqrt(3)
        row1_y_extent = [0, circle_diameter]
        row2_y_center = circle_radius + y_dist_rows
        row2_y_extent = [row2_y_center - circle_radius, row2_y_center + circle_radius]
        row3_y_center = circle_radius + 2 * y_dist_rows
        row3_y_extent = [row3_y_center - circle_radius, row3_y_center + circle_radius]

        # Calculate free space on the left side (x from 0 to 20)
        # This is a vertical gap between row 1 and row 3, as row 2 starts at x=20
        left_pocket_w = circle_radius
        left_pocket_h = row3_y_extent[0] - row1_y_extent[1]
        n_left = math.floor(left_pocket_w / square_side) * math.floor(left_pocket_h / square_side)
        print(f"Left-side pocket (size {left_pocket_w:.1f}x{left_pocket_h:.1f}) fits {n_left} squares.")
        N += n_left

        # Calculate free space on the right side (x from 120 to 140)
        # This space is split vertically by row 2
        right_pocket_w = sheet_w - (circle_radius * 2 + circle_radius * 2 * 2) # 140 - 120 = 20
        
        # Bottom-right pocket
        right_pocket_h_bottom = row2_y_extent[0]
        n_right_bottom = math.floor(right_pocket_w / square_side) * math.floor(right_pocket_h_bottom / square_side)
        print(f"Bottom-right pocket (size {right_pocket_w:.1f}x{right_pocket_h_bottom:.1f}) fits {n_right_bottom} squares.")
        N += n_right_bottom
        
        # Top-right pocket
        right_pocket_h_top = sheet_h - row2_y_extent[1]
        n_right_top = math.floor(right_pocket_w / square_side) * math.floor(right_pocket_h_top / square_side)
        print(f"Top-right pocket (size {right_pocket_w:.1f}x{right_pocket_h_top:.1f}) fits {n_right_top} squares.")
        N += n_right_top
    else:
        # Simplified case if M is not 9 (though we know it is)
        print("Optimal M is not 9, N calculation needs specific packing. Defaulting N=0.")
        
    print(f"The maximum number of squares (N) is {N}.")
    
    print("\n--- Step 5: Calculate Final Result ---")
    K = chars_per_square * N + chars_per_circle * M
    print(f"Based on the optimal values: N = {N}, M = {M}.")
    print(f"The final equation for the total characters is:")
    print(f"K = ({chars_per_square} * {N}) + ({chars_per_circle} * {M}) = {chars_per_square * N} + {chars_per_circle * M} = {K}")
    
    print("\nFinal Answer:")
    print(f"The ratio N:M:K is {N}:{M}:{K}")

calculate_max_artifacts()