import math

def solve_engraving_optimization():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """

    # Material dimensions in cm
    material_w = 140
    material_h = 110

    # Artifact properties
    square_dim = 10
    square_chars = 4
    
    # A circle of 20cm radius requires a 40x40cm square bounding box to be cut out
    circle_bounding_box_dim = 40
    circle_chars = 999

    # Step 1: Maximize the number of circles (M) first due to higher character value
    # We lay out the circle's 40x40cm bounding boxes in a grid
    circles_along_w = material_w // circle_bounding_box_dim
    circles_along_h = material_h // circle_bounding_box_dim
    M = circles_along_w * circles_along_h

    # Step 2: Calculate the area used by the circles
    used_area_w = circles_along_w * circle_bounding_box_dim
    used_area_h = circles_along_h * circle_bounding_box_dim

    # Step 3: Calculate the dimensions of the two leftover rectangular areas
    # Area 1: The strip remaining along the longer side of the material
    rem_area1_w = material_w - used_area_w
    rem_area1_h = material_h
    
    # Area 2: The strip remaining along the shorter side (on top of the used area)
    rem_area2_w = used_area_w
    rem_area2_h = material_h - used_area_h

    # Step 4: Calculate how many squares (N) can fit in the leftover areas
    squares_in_area1 = (rem_area1_w // square_dim) * (rem_area1_h // square_dim)
    squares_in_area2 = (rem_area2_w // square_dim) * (rem_area2_h // square_dim)
    N = squares_in_area1 + squares_in_area2
    
    # Step 5: Calculate the maximum total characters (K)
    K = (square_chars * N) + (circle_chars * M)

    # Print the results as requested
    print(f"Optimal number of squares (N): {N}")
    print(f"Optimal number of circles (M): {M}")
    print("\nTo find the maximum number of characters (K), we use the equation:")
    print(f"K = (Characters from Squares) + (Characters from Circles)")
    print(f"K = ({square_chars} * N) + ({circle_chars} * M)")
    print(f"K = ({square_chars} * {N}) + ({circle_chars} * {M})")
    print(f"K = {square_chars * N} + {circle_chars * M}")
    print(f"K = {K}")

    print(f"\nThe final answer in N:M:K format is: {N}:{M}:{K}")

solve_engraving_optimization()