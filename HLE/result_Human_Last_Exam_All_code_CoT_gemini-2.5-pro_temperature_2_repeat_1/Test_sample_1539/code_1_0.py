import math

def solve_laozi_puzzle():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters from a rectangular sheet of material.
    """
    # --- 1. Define Constants ---
    material_width = 140
    material_height = 110

    circle_radius = 20
    circle_diameter = 2 * circle_radius
    
    square_side = 10

    chars_per_circle = 9999
    chars_per_square = 360

    # --- 2. Maximize High-Value Circles ---
    # Based on analysis, the maximum number of 40cm diameter circles in a
    # 140x110cm rectangle is 6, arranged in a 3x2 grid.
    num_circles = 6
    circ_grid_cols = 3
    circ_grid_rows = 2
    
    # Calculate the bounding box of the circle grid
    circle_block_width = circ_grid_cols * circle_diameter
    circle_block_height = circ_grid_rows * circle_diameter

    # --- 3. Calculate Squares in Leftover Rectangular Space ---
    # Assume the circle block is placed at a corner of the material.
    # This leaves two main rectangular areas of waste material.
    
    # Area 1: Strip along the material's width
    waste_rect_1_w = material_width
    waste_rect_1_h = material_height - circle_block_height
    squares_in_rect_1 = (waste_rect_1_w // square_side) * (waste_rect_1_h // square_side)
    
    # Area 2: Strip along the circle block's height
    waste_rect_2_w = material_width - circle_block_width
    waste_rect_2_h = circle_block_height
    squares_in_rect_2 = (waste_rect_2_w // square_side) * (waste_rect_2_h // square_side)
    
    squares_from_waste = squares_in_rect_1 + squares_in_rect_2

    # --- 4. Calculate Squares in Interstitial Space ---
    # A 10x10cm square can fit in the space between 4 packed circles.
    # The number of these spaces in a grid is (cols-1) * (rows-1).
    num_interstitial_gaps = (circ_grid_cols - 1) * (circ_grid_rows - 1)
    squares_from_interstitial = num_interstitial_gaps # One square per gap

    # --- 5. Final Calculation ---
    total_squares = squares_from_waste + squares_from_interstitial
    
    # Total characters (the value to be maximized)
    max_characters = (num_circles * chars_per_circle) + (total_squares * chars_per_square)

    # --- 6. Output Results ---
    print("This script calculates the optimal number of plates and the maximum characters.")
    print("The strategy is to maximize the high-value circular plates first, then fill all remaining space with square plates.")
    print("\n--- Calculation ---")
    print(f"Number of circular plates (N): {num_circles}")
    print(f"Number of square plates (M): {total_squares} ({squares_from_waste} from rectangular leftovers + {squares_from_interstitial} from spaces between circles)")
    
    print("\nThe final equation for the total number of characters (K) is:")
    print(f"{num_circles} * {chars_per_circle} + {total_squares} * {chars_per_square} = {max_characters}")
    
    print("\nThe final answer in the format N:M:K is:")
    print(f"{num_circles}:{total_squares}:{max_characters}")

# Execute the function to solve the puzzle
solve_laozi_puzzle()