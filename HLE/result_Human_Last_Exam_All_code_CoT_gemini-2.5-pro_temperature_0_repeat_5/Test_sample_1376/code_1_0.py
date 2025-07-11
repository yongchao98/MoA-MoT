def solve_engraving_problem():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # --- Define constants based on the problem description ---

    # Material dimensions in cm
    material_width = 140
    material_height = 110

    # Artifact dimensions in cm
    square_side = 10
    # A circle of 20cm radius fits in a square of 40cm side
    circle_bounding_box_side = 40

    # Number of characters per artifact
    chars_per_square = 4
    chars_per_circle = 999

    # --- Step 1: Maximize the number of circles (M) ---
    # The value per area for circles is much higher, so we prioritize them.
    
    # Fit circles along the width and height of the material
    num_circles_along_width = material_width // circle_bounding_box_side
    num_circles_along_height = material_height // circle_bounding_box_side
    
    M = num_circles_along_width * num_circles_along_height

    # --- Step 2: Calculate the number of squares (N) in the remaining area ---
    
    # Calculate total area and the area used by circles
    total_area = material_width * material_height
    area_per_circle_bbox = circle_bounding_box_side ** 2
    area_used_by_circles = M * area_per_circle_bbox
    
    # The remaining area will be used for squares
    remaining_area = total_area - area_used_by_circles
    
    # Calculate how many squares can fit in the remaining area
    area_per_square = square_side ** 2
    N = remaining_area // area_per_square

    # --- Step 3: Calculate the maximum total characters (K) ---
    K = (chars_per_square * N) + (chars_per_circle * M)

    # --- Step 4: Print the results ---
    print("To maximize the number of engraved characters:")
    print(f"Number of squares (N): {N}")
    print(f"Number of circles (M): {M}")
    print("\nThe maximum number of characters (K) is calculated as:")
    print(f"K = ({chars_per_square} * {N}) + ({chars_per_circle} * {M}) = {K}")
    
    print("\nThe final answer in the format N:M:K is:")
    print(f"{N}:{M}:{K}")

solve_engraving_problem()
<<<58:6:6226>>>