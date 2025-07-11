import math

def solve_laozi_engraving():
    """
    Calculates the optimal number of circular and squared plates to maximize
    the number of engraved characters.
    """
    # Step 1: Define values and dimensions
    sheet_w = 140
    sheet_h = 110
    
    circle_radius = 20
    circle_diameter = 2 * circle_radius
    circle_value = 9999
    
    square_side = 10
    square_value = 360

    # Step 2 & 3: Strategy is to maximize circles first.
    # A simple grid packing (140/40 x 110/40 = 3x2=6) is not optimal.
    # A staggered 3-2-3 packing pattern fits 8 circles.
    # Let's verify the dimensions of this packing:
    # Width: 3 circles in a row with centers at x=20, 60, 100.
    # Bounding box width = 100 (center of last circle) + 20 (radius) = 120 cm.
    # Height: 3 rows are staggered.
    # Vertical distance between staggered rows = sqrt(diameter^2 - radius^2)
    # This is incorrect. The horizontal distance between centers is diameter, but here they are packed tighter.
    # If centers are (20,y1), (60,y1), (100,y1) and staggered row is (40, y2), (80, y2),
    # the distance between centers (20,y1) and (40,y2) is 40.
    # vert_dist = sqrt(40^2 - 20^2) = sqrt(1200) ~= 34.64 cm.
    # Total height = radius_1st_row + vert_dist + vert_dist + radius_3rd_row
    # height = 20 + 34.64 + 34.64 + 20 = 109.28 cm.
    # This 120cm x 109.28cm packing fits within the 140x110cm sheet.
    
    num_circles = 8
    
    # We allocate a 120x110 cm area for the 8 circles.
    circles_area_w = 120
    circles_area_h = 110

    # Step 4: Calculate remaining material and the number of squares.
    remaining_w = sheet_w - circles_area_w
    remaining_h = sheet_h
    
    # The remaining area is a rectangle of 20x110 cm.
    num_squares_in_w = math.floor(remaining_w / square_side)
    num_squares_in_h = math.floor(remaining_h / square_side)
    num_squares = num_squares_in_w * num_squares_in_h

    # Step 5: Final calculation
    N = num_circles
    M = num_squares
    K = (N * circle_value) + (M * square_value)
    
    # Print the equation as requested
    print(f"{N} * {circle_value} + {M} * {square_value} = {K}")
    
    # Print the final answer in the N:M:K format
    final_answer = f"{N}:{M}:{K}"
    print(f"\nThe optimal production is {N} circular plates and {M} squared plates.")
    print(f"The maximum number of characters is {K}.")
    print(f"\nFinal Answer Format: {final_answer}")
    
    return final_answer

# Execute the function
final_answer_string = solve_laozi_engraving()

# The final line formatted for the platform
print(f"\n<<<"+final_answer_string+">>>")
