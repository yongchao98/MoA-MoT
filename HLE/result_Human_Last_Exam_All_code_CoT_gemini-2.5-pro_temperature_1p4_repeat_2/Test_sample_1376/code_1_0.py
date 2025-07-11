import math

# Step 1: Define constants based on the problem description.
MATERIAL_WIDTH = 140
MATERIAL_HEIGHT = 110

SQUARE_SIDE = 10
CHARS_PER_SQUARE = 4

CIRCLE_RADIUS = 20
# For packing, a circle is treated as its bounding square (diameter x diameter).
CIRCLE_DIAMETER = CIRCLE_RADIUS * 2
CHARS_PER_CIRCLE = 999

def solve_engraving_maximization():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # Step 2: Determine the maximum number of circles (M).
    # The strategy is to prioritize circles as they yield more characters per area.
    # We check how many 40x40cm bounding boxes for circles fit into the material.
    
    # An 80x120cm block (2x3 circles) would not fit, as 120 > 110.
    # Therefore, the only possible arrangement for maximum circles is 3x2.
    num_circles_along_width = MATERIAL_WIDTH // CIRCLE_DIAMETER  # 140 // 40 = 3
    num_circles_along_height = MATERIAL_HEIGHT // CIRCLE_DIAMETER # 110 // 40 = 2
    
    M = num_circles_along_width * num_circles_along_height

    # Calculate the dimensions of the rectangular block occupied by the packed circles.
    circles_block_width = num_circles_along_width * CIRCLE_DIAMETER
    circles_block_height = num_circles_along_height * CIRCLE_DIAMETER

    # Step 3: Calculate the maximum number of squares (N) in the remaining area.
    # The remaining area is an L-shape if we place the circle block in a corner.
    # We can calculate the squares in this L-shape by splitting it into two rectangles.
    
    # Rectangle 1 of the L-shape
    rem_rect1_w = MATERIAL_WIDTH - circles_block_width # 140 - 120 = 20
    rem_rect1_h = MATERIAL_HEIGHT                   # 110
    squares_in_rect1 = (rem_rect1_w // SQUARE_SIDE) * (rem_rect1_h // SQUARE_SIDE) # (2*11)=22

    # Rectangle 2 of the L-shape
    rem_rect2_w = circles_block_width                # 120
    rem_rect2_h = MATERIAL_HEIGHT - circles_block_height # 110 - 80 = 30
    squares_in_rect2 = (rem_rect2_w // SQUARE_SIDE) * (rem_rect2_h // SQUARE_SIDE) # (12*3)=36

    N = squares_in_rect1 + squares_in_rect2

    # Step 4: Calculate the total number of engraved characters (K).
    K = (N * CHARS_PER_SQUARE) + (M * CHARS_PER_CIRCLE)

    # Step 5: Print the detailed results and the final answer.
    print(f"To maximize engraved characters, the optimal numbers of artifacts are:")
    print(f"Number of Squares (N): {N}")
    print(f"Number of Circles (M): {M}")
    print("\nThe total number of characters (K) is calculated as:")
    print(f"K = (N squares * {CHARS_PER_SQUARE} chars/square) + (M circles * {CHARS_PER_CIRCLE} chars/circle)")
    print(f"K = ({N} * {CHARS_PER_SQUARE}) + ({M} * {CHARS_PER_CIRCLE})")
    print(f"K = {N * CHARS_PER_SQUARE} + {M * CHARS_PER_CIRCLE}")
    print(f"K = {K}")
    
    print("\nThe final answer in the format N:M:K is:")
    print(f"{N}:{M}:{K}")

if __name__ == '__main__':
    solve_engraving_maximization()