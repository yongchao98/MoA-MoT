import math

def solve_emperor_tomb_problem():
    """
    Calculates the optimal number of squares and circles to maximize the number of engraved characters.
    """
    # --- Step 1: Define constants from the problem description ---
    material_w = 140
    material_h = 110

    square_size = 10
    square_chars = 4

    circle_radius = 20
    circle_diameter = circle_radius * 2
    circle_chars = 999

    # --- Step 2: Maximize the number of circles (M) ---
    # We use a greedy approach, prioritizing circles due to their high character count.
    # We find the maximum number of circles that fit in a grid layout.
    
    # Check orientation 1: 140x110
    num_circles_w1 = math.floor(material_w / circle_diameter)
    num_circles_h1 = math.floor(material_h / circle_diameter)
    M1 = num_circles_w1 * num_circles_h1
    
    # Check orientation 2: 110x140 (by swapping width and height)
    num_circles_w2 = math.floor(material_h / circle_diameter)
    num_circles_h2 = math.floor(material_w / circle_diameter)
    M2 = num_circles_w2 * num_circles_h2
    
    # The max number of circles is the same in both orientations
    M = max(M1, M2)
    # The layout for the circles will be 3x2.
    num_circles_w = 3
    num_circles_h = 2

    print(f"The emperor's workers can cut a maximum of {M} circles.")

    # --- Step 3: Calculate the number of squares (N) in the remaining area ---
    # The block of 6 circles will occupy a (3*40)x(2*40) = 120x80 cm area.
    circles_block_w = num_circles_w * circle_diameter
    circles_block_h = num_circles_h * circle_diameter
    
    # By placing this block in the corner of the 140x110 cm material,
    # we are left with two rectangular areas.
    
    # Area 1: A strip along the bottom
    rem_area1_w = material_w
    rem_area1_h = material_h - circles_block_h # 110 - 80 = 30
    
    # Area 2: A strip down the side
    rem_area2_w = material_w - circles_block_w # 140 - 120 = 20
    rem_area2_h = circles_block_h # This strip is next to the 80cm side of the circle block

    # Calculate how many 10x10 squares fit in each remaining rectangle
    squares_in_area1 = math.floor(rem_area1_w / square_size) * math.floor(rem_area1_h / square_size)
    squares_in_area2 = math.floor(rem_area2_w / square_size) * math.floor(rem_area2_h / square_size)
    
    N = squares_in_area1 + squares_in_area2
    
    print(f"The remaining space can be filled with {N} squares.")

    # --- Step 4: Calculate the maximum total characters (K) ---
    # K = (characters per square * N) + (characters per circle * M)
    K = (square_chars * N) + (circle_chars * M)

    print("\nCalculating the maximum number of characters (K):")
    # Output the final equation with all numbers
    print(f"K = {square_chars} * {N} + {circle_chars} * {M} = {square_chars * N} + {circle_chars * M} = {K}")

    # --- Step 5: Format the final answer ---
    final_answer = f"{N}:{M}:{K}"
    print(f"\nThe final answer in the format N:M:K is: {final_answer}")
    
    return final_answer

# Execute the function and capture the final answer for the required format
final_answer_value = solve_emperor_tomb_problem()
# The final response should end with this format
# print(f"<<<{final_answer_value}>>>") # This is for the system, not the user output

if __name__ == '__main__':
    # This block is for direct execution and is not part of the final assistant response.
    # The required output format is handled by the thought process.
    pass
<<<58:6:6226>>>