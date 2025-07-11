import math

# A simple class for 2D points to make calculations cleaner.
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f})"

def main():
    """
    Solves the visibility puzzle by modeling the room and performing line-of-sight checks.
    """
    print("Analyzing Joe's visibility in the room.")
    print("A 2D coordinate system is established with the southwest corner at (0, 0).")
    print("The room is a 12x12 foot square.\n")

    # Define Joe's position
    joe_x_min = 4.5
    joe_x_max = 7.5
    joe_y_max = 0.25
    joe_eye_height = 5.0
    print(f"Joe's viewing position is a rectangle from x={joe_x_min} to {joe_x_max} and y=0 to {joe_y_max}.")
    print(f"Joe's eye height is {joe_eye_height} feet.\n")

    # --- Define Obstacles based on the room description ---
    # We assume the door (3ft wide) is opened into the room, perpendicular to the south wall.
    # We assume the hinge is on the left (at x=4.5), so the door itself forms a vertical barrier.
    door_hinge = Point(4.5, 0)
    door_edge = Point(4.5, 3)
    bookshelf_corner = Point(1, 4)
    wardrobe_body_sw = Point(9.5, 4)
    wardrobe_body_ne = Point(12, 8)
    wardrobe_door1_corner = Point(7.5, 4)
    wardrobe_door2_corner = Point(7.5, 8)
    wardrobe_front_corner_top = Point(9.5, 8)
    
    print("Obstacle setup:")
    print(f" - Door: A 3-foot wall segment from {door_hinge} to {door_edge}.")
    print(f" - Bookshelf: An obstacle from (0,0) to {bookshelf_corner}.")
    print(f" - Wardrobe: A solid body from {wardrobe_body_sw} to {wardrobe_body_ne}.")
    print(f" - Wardrobe Doors: Lines from {wardrobe_body_sw} to {wardrobe_door1_corner} and from {wardrobe_front_corner_top} to {wardrobe_door2_corner}.\n")
    
    visible_balls = []

    # --- Analysis for each ball ---

    print("--- Red Ball (SE corner) ---")
    ball_red = Point(11.75, 0.25)
    print(f"The red ball is on the floor at {ball_red}.")
    # Joe can look along the south wall. Check from a valid viewing point.
    viewer = Point(7.5, 0.25)
    print(f"From viewing point {viewer}, the line of sight to the ball is along the line y={viewer.y}.")
    print("This path is not blocked by the bookshelf (x<1), wardrobe (y>4), or wardrobe doors (y>4).")
    print(f"The door at x={door_hinge.x} is behind this particular viewpoint.")
    print("Result: Red ball is VISIBLE.\n")
    visible_balls.append("Red")

    print("--- Green Ball (SW corner, on top of bookshelf) ---")
    ball_green_xy = Point(0.25, 0.25)
    ball_green_z = 7.25
    bookshelf_height = 7.0
    bookshelf_front_x = 1.0
    print(f"The green ball is at {ball_green_xy} on a {bookshelf_height} ft tall shelf. Ball center height is {ball_green_z} ft.")
    
    print("Step 1: Check if Joe can see *over* the 7-foot bookshelf.")
    # Condition from analysis: vx > 7 for the line of sight to be high enough.
    vx_for_3d_clearance = 7.0
    print(f"3D calculation shows Joe must stand at an x-position vx > {vx_for_3d_clearance:.2f} to see over the shelf.")
    print(f"This is possible, as Joe's viewing area extends to x=7.5.")
    
    print("Step 2: Check for 2D obstructions from that position, primarily the door at x=4.5.")
    viewer_vx = 7.1
    ball_bx = 0.25
    viewer_vy = 0.25
    ball_by = 0.25
    y_intersect_at_door = viewer_vy + ((ball_by - viewer_vy) / (ball_bx - viewer_vx)) * (door_hinge.x - viewer_vx)
    print(f"A viewer at ({viewer_vx}, {viewer_vy}) looking at the ball at ({ball_bx}, {ball_by}) has a line of sight that crosses the door plane (x={door_hinge.x}).")
    print(f"The intersection y-coordinate is: {y_intersect_at_door:.2f}.")
    print(f"The physical door occupies the space from y=0 to y={door_edge.y}.")
    print(f"Since y={y_intersect_at_door:.2f} is within the door's range [0, {door_edge.y}], the view is blocked.")
    print("Result: Green ball is NOT VISIBLE.\n")

    print("--- Yellow Ball (NW corner) ---")
    ball_yellow = Point(0.25, 11.75)
    print(f"The yellow ball is on the floor at {ball_yellow}.")
    # To see the yellow ball, Joe needs to stand on the right side of the doorway.
    viewer = Point(7.5, 0.25)
    print(f"The best viewing position is the far right, {viewer}.")
    m = (ball_yellow.y - viewer.y) / (ball_yellow.x - viewer.x)
    
    y_at_door = viewer.y + m * (door_hinge.x - viewer.x)
    print(f"The line of sight intersects the door plane (x={door_hinge.x}) at y = {y_at_door:.2f}.")
    print(f"This is above the door's top edge (y={door_edge.y}), so the door doesn't block the view.")
    
    y_at_shelf = viewer.y + m * (bookshelf_corner.x - viewer.x)
    print(f"The line of sight passes the bookshelf's front plane (x={bookshelf_corner.x}) at y = {y_at_shelf:.2f}.")
    print(f"This is beyond the bookshelf's north edge (y={bookshelf_corner.y}), so the bookshelf doesn't block the view.")
    print("Result: Yellow ball is VISIBLE.\n")
    visible_balls.append("Yellow")

    print("--- Purple Ball (in SE corner of wardrobe) ---")
    ball_purple = Point(11.75, 4.25)
    print(f"The purple ball is inside the wardrobe, at {ball_purple}.")
    print(f"The wardrobe front is on the line x={wardrobe_body_sw.x}, with an opening from y={wardrobe_body_sw.y} to y={wardrobe_body_ne.y}.")
    # To see the purple ball, the best chance is from the far left of the doorway.
    viewer = Point(4.5, 0.25)
    print(f"Checking from the best viewpoint, approximately {viewer}.")
    m = (ball_purple.y - viewer.y) / (ball_purple.x - viewer.x)
    y_at_wardrobe = viewer.y + m * (wardrobe_body_sw.x - viewer.x)
    print(f"The line of sight from {viewer} hits the wardrobe front (x={wardrobe_body_sw.x}) at y = {y_at_wardrobe:.2f}.")
    print(f"Since {y_at_wardrobe:.2f} is below the opening's start at y={wardrobe_body_sw.y}, the view is blocked.")
    print("Result: Purple ball is NOT VISIBLE.\n")

    print("--- Blue Ball (NE corner) ---")
    ball_blue = Point(11.75, 11.75)
    print(f"The blue ball is at {ball_blue}.")
    print("The main potential obstruction is the corner of the wardrobe and its open door, near point (9.5, 8).")
    print("For a clear view, the line of sight must pass this corner. This leads to an inequality for the viewer's position (vx).")
    print("   slope_to_ball > slope_to_obstacle_corner")
    print("   11.5 / (11.75 - vx) > 7.75 / (9.5 - vx)")
    
    vx_max_for_blue = 18.1875 / 3.75
    print(f"Solving the inequality for vx gives: vx < {vx_max_for_blue:.2f}")
    print(f"Joe can stand at any x-position from {joe_x_min} to {joe_x_max}.")
    print(f"Since there is an overlap (e.g., from x={joe_x_min} to x={vx_max_for_blue:.2f}), Joe can find a spot to see the ball.")
    print("Result: Blue ball is VISIBLE.\n")
    visible_balls.append("Blue")

    print("------ Summary ------")
    print("Based on the geometric analysis, the balls Joe can see are:")
    # Sort for consistent output
    final_visible_list = sorted(visible_balls)
    for ball in final_visible_list:
        print(f"  - {ball}")

    # Final answer format
    final_answer = ", ".join(final_visible_list)
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    main()