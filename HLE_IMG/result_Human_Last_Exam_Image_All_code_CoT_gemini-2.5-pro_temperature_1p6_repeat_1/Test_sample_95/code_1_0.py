import math

def solve_puzzle():
    """
    Solves the AGI QA game to find the center of the right-most white circle.
    """

    # Step 1 & 2: Determine the radii of the circles
    # Diameter of yellow circle is 1 cm
    d_y = 1.0
    # Radius of yellow circle
    r_y = d_y / 2
    
    print(f"Given the diameter of a yellow circle is {d_y} cm, its radius (r_y) is {r_y} cm.")

    # A key geometric configuration is a small circle (yellow) tangent to a line (boundary)
    # and two larger tangent circles (white). This leads to the relationship R = 4r.
    # Let's derive this:
    # Let white circle radius be R and yellow circle radius be r.
    # Assume white circles are tangent to the boundary.
    # Place two white circle centers at (R, R) and (3R, R) and the boundary at y=0. This isn't the config.
    # Correct config: White circles centers at (x1, R), (x2, R) where x2-x1=2R. Boundary at y=2R.
    # Yellow circle center is at ((x1+x2)/2, 2R-r).
    # Distance from yellow center to white center must be R+r.
    # ( (x2-x1)/2 )^2 + ( (2R-r) - R )^2 = (R+r)^2
    # R^2 + (R-r)^2 = (R+r)^2
    # R^2 + R^2 - 2Rr + r^2 = R^2 + 2Rr + r^2
    # R^2 = 4Rr
    # R = 4r
    
    R_w = 4 * r_y
    print(f"From the geometry of the circles and the boundary, we deduce that the radius of a white circle (R_w) is 4 times the radius of a yellow circle (r_y).")
    print(f"R_w = 4 * {r_y} cm = {R_w} cm.")
    # This aligns with the 0.5 cm grid rule.

    # Step 3: Determine overall layout and dimensions
    # The layout is a 4x3 grid of white circles.
    num_columns = 4
    num_rows = 3
    
    # The total width is determined by 4 circles in a row, tangent to each other and the boundaries.
    # Width = R_w (left edge to center1) + 3 * (2*R_w) (center1 to center4) + R_w (center4 to right edge) = 8*R_w
    W = 8 * R_w
    # The total height is determined by 3 circles in a column.
    # Height = R_w (bottom edge to center1) + 2 * (2*R_w) (center1 to center3) + R_w (center3 to top edge) = 6*R_w
    H = 6 * R_w
    
    print(f"The layout is a {num_columns}x{num_rows} grid of white circles with radius {R_w} cm.")
    print(f"The total width of the area is 8 * R_w = 8 * {R_w} = {W} cm.")
    print(f"The total height of the area is 6 * R_w = 6 * {R_w} = {H} cm.")
    
    # Step 4 & 5: Establish coordinate system and find coordinates
    # Let origin (0,0) be the bottom-left corner.
    # x-coordinates of the centers of the 4 columns:
    x_coords = [R_w + i * (2 * R_w) for i in range(num_columns)]
    # y-coordinates of the centers of the 3 rows:
    y_coords = [R_w + i * (2 * R_w) for i in range(num_rows)]
    
    # The right-most column of circles is the last one.
    rightmost_x = x_coords[-1]
    
    # The question asks for "the" right-most white circle (singular).
    # This suggests picking one from the right-most column. The middle one is the most logical choice.
    # The y-coordinate for the middle row is y_coords[1].
    middle_row_y = y_coords[1]
    
    # Final coordinates for the center of the target circle
    target_x = rightmost_x
    target_y = middle_row_y
    
    print(f"\nPlacing the origin (0,0) at the bottom-left corner:")
    print(f"The x-coordinate of the center of the right-most column of circles is {target_x} cm.")
    print(f"The y-coordinate of the center of the middle row of circles is {target_y} cm.")
    print("\nThe question asks for the center of 'the right-most white circle'. The most distinct circle in the right-most column is the one in the middle row.")
    print(f"Therefore, the center of the right-most white circle is at x = {target_x}, y = {target_y}.")
    
    # Present final answer in requested format
    final_answer = f"{int(target_x)}:{int(target_y)}"
    print(f"\nThe final answer is...")

solve_puzzle()
# The final answer is calculated based on the derived geometry.
print("14:6")