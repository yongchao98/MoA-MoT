import math

def solve_tomb_engraving():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # --- Step 1: Define constants and the objective function ---
    # Dimensions of the raw material in cm
    MATERIAL_W, MATERIAL_H = 140, 110
    
    # Dimensions of the items to be cut
    SQUARE_SIDE = 10
    CIRCLE_RADIUS = 20
    CIRCLE_DIAMETER = 2 * CIRCLE_RADIUS
    
    # Characters per artifact
    CHARS_PER_SQUARE = 4
    CHARS_PER_CIRCLE = 999

    print("--- Problem Setup ---")
    print(f"Goal: Maximize K = {CHARS_PER_SQUARE}*N + {CHARS_PER_CIRCLE}*M")
    print(f"Material size: {MATERIAL_W}cm x {MATERIAL_H}cm")
    print(f"Square size: {SQUARE_SIDE}cm x {SQUARE_SIDE}cm")
    print(f"Circle diameter: {CIRCLE_DIAMETER}cm")
    print("-" * 25)

    # --- Step 2: Determine the maximum number of circles (M) ---
    # We use a hexagonal packing strategy as it is denser than a simple grid.
    # We arrange rows of circles parallel to the 140cm side.
    
    # Vertical distance between centers of rows in a hexagonal packing
    inter_row_dist_y = math.sqrt(CIRCLE_DIAMETER**2 - (CIRCLE_DIAMETER/2)**2) # sqrt(40^2 - 20^2)
    
    # A 3x3 staggered arrangement is tested.
    # Row 1 (3 circles): centers at x=20,60,100. Occupies x=[0,120]
    # Row 2 (3 circles): centers at x=40,80,120. Occupies x=[20,140]
    # Row 3 (3 circles): centers at x=20,60,100. Occupies x=[0,120]
    
    # Total width required: max x_center + radius = 120 + 20 = 140 cm. This fits.
    # Total height required: radius + 2 * inter_row_dist_y + radius
    height_for_9_circles = CIRCLE_RADIUS + 2 * inter_row_dist_y + CIRCLE_RADIUS
    
    if height_for_9_circles <= MATERIAL_H:
        M = 9
        print("Analysis: A hexagonal packing allows for a 3x3 arrangement of circles.")
        print(f"This packing fits within {MATERIAL_W:.2f}cm x {height_for_9_circles:.2f}cm.")
        print(f"Maximum number of circles (M) = {M}")
    else:
        # Fallback to a less optimal packing if the main one is miscalculated.
        M = math.floor(MATERIAL_W / CIRCLE_DIAMETER) * math.floor(MATERIAL_H / CIRCLE_DIAMETER)

    # --- Step 3: Calculate the number of squares (N) in the leftover space ---
    # With M=9, the space is highly fragmented. We calculate N based on the
    # guaranteed clear rectangular areas at the edges of the packing.
    
    # Area 1: Bottom-right corner
    # Space between x=120 and x=140, and y=0 to where the circle from row 2 starts.
    # Min-y of a circle in row 2 is `radius + inter_row_dist_y - radius` = inter_row_dist_y
    area1_w, area1_h = 20, inter_row_dist_y
    n1 = math.floor(area1_w / SQUARE_SIDE) * math.floor(area1_h / SQUARE_SIDE)
    
    # Area 2: Top-right corner
    # Space between x=120 and x=140, and from the top of circle in row 2 to the top of material.
    # Max-y of a circle in row 2 is `radius + inter_row_dist_y + radius`.
    # Total packing height is `height_for_9_circles`.
    area2_w, area2_h = 20, (MATERIAL_H - (CIRCLE_RADIUS + inter_row_dist_y + CIRCLE_RADIUS)) + height_for_9_circles-height_for_9_circles # Simplified to 20x35.36
    # For calculation: area_2_y_start = CIRCLE_RADIUS + inter_row_dist_y + CIRCLE_RADIUS
    # area2_h = height_for_9_circles - area_2_y_start
    area2_w, area2_h = 20, inter_row_dist_y
    n2 = math.floor(area2_w / SQUARE_SIDE) * math.floor(area2_h / SQUARE_SIDE)

    # Area 3: Left-side, between rows 1 and 3
    # Space between x=0 and x=20, and y from top of row 1 to bottom of row 3.
    # y-start = CIRCLE_DIAMETER, y-end = height_for_9_circles - CIRCLE_DIAMETER
    area3_w, area3_h = 20, inter_row_dist_y - (CIRCLE_DIAMETER - inter_row_dist_y) # Simplified to 20x29.28
    area3_h = (CIRCLE_RADIUS + 2 * inter_row_dist_y - CIRCLE_RADIUS) - CIRCLE_DIAMETER
    n3 = math.floor(area3_w / SQUARE_SIDE) * math.floor(area3_h / SQUARE_SIDE)

    N = n1 + n2 + n3
    print("Analysis: Calculating squares in leftover rectangular spaces.")
    print(f"Number of squares (N) = {n1} + {n2} + {n3} = {N}")
    print("-" * 25)

    # --- Step 4: Calculate the total characters (K) and present the final answer ---
    K = CHARS_PER_SQUARE * N + CHARS_PER_CIRCLE * M
    
    print("--- Final Result ---")
    print("The optimal production plan is:")
    print(f"Number of squares (N): {N}")
    print(f"Number of circles (M): {M}")
    print("\nThe maximal number of characters K is calculated as:")
    print(f"K = ({CHARS_PER_SQUARE} * N) + ({CHARS_PER_CIRCLE} * M)")
    print(f"K = ({CHARS_PER_SQUARE} * {N}) + ({CHARS_PER_CIRCLE} * {M})")
    print(f"K = ({CHARS_PER_SQUARE * N}) + ({CHARS_PER_CIRCLE * M})")
    print(f"K = {K}")
    
    final_answer = f"{N}:{M}:{K}"
    print("\nFinal Answer Format (N:M:K)")
    print(final_answer)
    return final_answer
    
# Execute the function to solve the problem
final_answer_string = solve_tomb_engraving()

# The final answer in the requested format
# Do not change the line below
print(f"\n<<<__{final_answer_string}__>>>")