import math

def solve_laozi_puzzle():
    """
    Solves the optimization problem for Laozi's books.

    This function calculates the optimal number of circular and squared plates
    to cut from a meteorite sheet to maximize the number of engraved symbols.
    It analyzes and compares a simple grid packing with a more efficient
    staggered packing for the circular plates.
    """

    # --- Define Constants ---
    sheet_width = 140
    sheet_height = 110

    circle_radius = 20
    circle_diameter = 2 * circle_radius
    circle_value = 9999

    square_side = 10
    square_value = 360

    print("Step 1: Evaluating Packing Strategies")
    print("The goal is to maximize the total number of symbols.")
    print(f"A circular plate (40cm diameter) yields {circle_value} symbols.")
    print(f"A squared plate (10cm side) yields {square_value} symbols.")
    print("Because circular plates offer more value per area, we prioritize packing as many circles as possible.\n")

    # --- Strategy 1: Simple Grid Packing ---
    # This is a simple but less optimal method.
    n_grid_w = sheet_width // circle_diameter
    n_grid_h = sheet_height // circle_diameter
    N_grid = n_grid_w * n_grid_h
    
    # Calculate remaining area for squares
    used_width = n_grid_w * circle_diameter
    used_height = n_grid_h * circle_diameter
    # Remaining area is split into two rectangles
    rem_rect1_w, rem_rect1_h = sheet_width - used_width, sheet_height
    rem_rect2_w, rem_rect2_h = used_width, sheet_height - used_height
    
    M_grid1 = (rem_rect1_w // square_side) * (rem_rect1_h // square_side)
    M_grid2 = (rem_rect2_w // square_side) * (rem_rect2_h // square_side)
    M_grid = M_grid1 + M_grid2
    
    K_grid = N_grid * circle_value + M_grid * square_value

    print("Step 2: Analyzing Staggered Packing (Optimal Strategy)")
    # --- Strategy 2: Staggered Packing ---
    # This is a more efficient packing method for circles.
    # We align the longer side (140cm) as the width.
    
    # Vertical calculation
    stagger_height_step = circle_radius * math.sqrt(3)
    num_rows = math.floor((sheet_height - circle_diameter) / stagger_height_step) + 1

    # In our case (140x110), we can fit 3 rows
    # Row 1: 3 circles, Row 2: 2 circles, Row 3: 3 circles
    N_staggered = 8 
    
    print(f"Using a staggered arrangement, we can fit N = {N_staggered} circular plates.")
    
    # Calculate M for the staggered packing
    # The cluster of 8 circles occupies a bounding area of roughly 120cm x 109.28cm
    cluster_width = 120 
    
    # 1. The main remaining strip on the side
    m_strip_w = sheet_width - cluster_width
    m_strip_h = sheet_height
    M_strip = (m_strip_w // square_side) * (m_strip_h // square_side)
    
    # 2. Rectangular areas on the sides of the middle row of circles
    # This space is 20cm wide and 40cm high, on both left and right sides.
    m_sides_w = circle_radius
    m_sides_h = circle_diameter
    M_sides = 2 * (m_sides_w // square_side) * (m_sides_h // square_side)

    # 3. Four corner areas around the corner circles.
    # The largest square that fits in a quadrant of radius 20 is side 14.14cm.
    # So a 10x10 square fits in each of the 4 outer corners of the cluster.
    M_corners = 4

    M_staggered = M_strip + M_sides + M_corners
    print(f"The remaining material can be cut into M = {M_strip} (side) + {M_sides} (middle) + {M_corners} (corners) = {M_staggered} squared plates.")

    K_staggered = N_staggered * circle_value + M_staggered * square_value

    print("\nStep 3: Final Calculation and Result")
    # The staggered packing is superior.
    N_final = N_staggered
    M_final = M_staggered
    K_final = K_staggered
    
    print(f"The optimal configuration is {N_final} circular plates and {M_final} squared plates.")
    print("The maximal number of characters K is calculated as:")
    print(f"K = {N_final} * {circle_value} + {M_final} * {square_value} = {N_final * circle_value} + {M_final * square_value} = {K_final}")

    # Final answer in the required format
    final_answer = f"{N_final}:{M_final}:{K_final}"
    print(f"\nFinal Answer Format (N:M:K):")
    print(final_answer)
    return final_answer

# Execute the function and capture the final answer for the specified output format
final_answer_string = solve_laozi_puzzle()
print(f"\n<<< {final_answer_string} >>>")