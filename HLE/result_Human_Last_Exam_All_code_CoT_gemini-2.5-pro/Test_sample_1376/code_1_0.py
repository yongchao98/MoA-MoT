import math

def solve_tomb_artifacts():
    """
    Calculates the optimal number of squares and circles to maximize the number of engraved characters.
    """
    # --- Define Constants ---
    material_w = 140
    material_h = 110

    square_size = 10
    circle_radius = 20
    circle_diameter = circle_radius * 2

    chars_per_square = 4
    chars_per_circle = 999

    # --- Strategy: Maximize the most valuable item (circles) first ---
    
    # We treat circles as their 40x40cm bounding boxes for packing
    
    # Calculate how many circles fit along each dimension of the material
    circles_along_w = math.floor(material_w / circle_diameter)
    circles_along_h = math.floor(material_h / circle_diameter)
    
    # The maximum number of circles (M) is their product
    M = circles_along_w * circles_along_h
    
    # --- Calculate the number of squares (N) from the remaining material ---
    
    # Calculate the dimensions of the area occupied by the grid of circles
    occupied_w = circles_along_w * circle_diameter
    occupied_h = circles_along_h * circle_diameter
    
    # The remaining area can be seen as two rectangular strips
    # Strip 1:
    strip1_w = material_w - occupied_w
    strip1_h = material_h
    squares_in_strip1 = math.floor(strip1_w / square_size) * math.floor(strip1_h / square_size)
    
    # Strip 2:
    strip2_w = occupied_w
    strip2_h = material_h - occupied_h
    squares_in_strip2 = math.floor(strip2_w / square_size) * math.floor(strip2_h / square_size)

    # Total number of squares (N)
    N = squares_in_strip1 + squares_in_strip2

    # --- Calculate the total number of characters (K) ---
    K = (chars_per_square * N) + (chars_per_circle * M)

    # --- Print the results ---
    print("Based on the optimization, the plan is as follows:")
    print(f"1. Cut a {circles_along_w}x{circles_along_h} grid of circles from the material.")
    print(f"   - Optimal number of circles (M): {M}")
    print("2. Cut squares from the two remaining rectangular strips of material.")
    print(f"   - Optimal number of squares (N): {N}")

    print("\nThis configuration maximizes the total number of characters.")
    print("\nThe final calculation for the maximal number of characters (K) is:")
    print(f"{chars_per_square} * {N} (from squares) + {chars_per_circle} * {M} (from circles) = {K}")

    print("\nThe final answer in the format N:M:K is:")
    print(f"{N}:{M}:{K}")


solve_tomb_artifacts()