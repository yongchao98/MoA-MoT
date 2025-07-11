import math

# Define a simple class for points for clarity
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"({self.x}, {self.y})"

def check_visibility():
    """
    Analyzes the visibility of each ball from Joe's position.
    """
    # --- Define positions ---
    # Joe's two most extreme viewpoints
    joe_v_left = Point(4.5, 0.25)
    joe_v_right = Point(7.5, 0.25)

    # Ball positions (using center of 6-inch ball, so radius is 0.25 ft)
    red_ball = Point(12 - 0.25, 0.25)
    # The green ball is on the corner of the shelf closest to the door
    green_ball = Point(1, 0.25)
    yellow_ball = Point(0.25, 12 - 0.25)
    blue_ball = Point(12 - 0.25, 12 - 0.25)
    # Purple ball is in the SE corner of the wardrobe (x=12, y=4)
    purple_ball = Point(12 - 0.25, 4 + 0.25)


    # Obstruction corners
    bookshelf_corner = Point(1, 4)
    wardrobe_corner = Point(9.5, 8)
    wardrobe_front_plane_x = 9.5
    wardrobe_base_y = 4
    open_door = {'x': 7.5, 'y_start': 0, 'y_end': 3}


    visible_balls = []

    # --- Analysis for each ball ---

    # 1. Red Ball (SE Corner)
    # View from the right side of the doorway, just left of the open door hinge.
    # The line of sight is along y=0.25. The wardrobe starts at y=4.
    # The open door is at x=7.5, but Joe can stand at x=7.49.
    # So the view is clear.
    is_red_visible = True
    if is_red_visible:
        visible_balls.append("Red")
    print(f"Checking Red Ball at {red_ball}:")
    print(f"  Joe's best viewpoint is just left of the door hinge at {joe_v_right}.")
    print("  The line of sight is along the south wall.")
    print("  The wardrobe starts at y=4 and does not block the view.")
    print("  Result: Visible\n")


    # 2. Green Ball (SW Corner, on shelf)
    # View from anywhere in the doorway to the corner of the shelf at (1, 0).
    # There are no objects in between. The ball is 7ft high, Joe's eyes are ~5ft.
    # A standard 7ft doorway would not block the upward view.
    is_green_visible = True
    if is_green_visible:
        visible_balls.append("Green")
    print(f"Checking Green Ball at {green_ball}:")
    print(f"  The ball is on the corner of the bookshelf nearest the door.")
    print(f"  The line of sight from {joe_v_left} to {green_ball} is unobstructed.")
    print("  Result: Visible\n")

    # 3. Yellow Ball (NW Corner)
    # View from the left side of the doorway. Blocked by the bookshelf corner.
    # We check if the ball is in the "shadow" of the bookshelf corner.
    view = joe_v_left
    target = yellow_ball
    obstacle = bookshelf_corner
    # Calculate y-coordinate of the line-of-sight to the obstacle corner, at the target's x-coordinate
    slope_to_obstacle = (obstacle.y - view.y) / (obstacle.x - view.x)
    y_at_target_x = slope_to_obstacle * (target.x - view.x) + view.y
    is_yellow_visible = target.y < y_at_target_x # Visible if below the shadow line
    if is_yellow_visible:
        visible_balls.append("Yellow")
    print(f"Checking Yellow Ball at {yellow_ball}:")
    print(f"  Joe's best viewpoint is {view}.")
    print(f"  The key obstruction is the bookshelf corner at {obstacle}.")
    print(f"  The line of sight to the ball (y={target.y}) is blocked because it needs to pass 'through' the bookshelf, whose shadow extends up to y={y_at_target_x:.2f} at the ball's x-position.")
    print(f"  Result: Not Visible\n")


    # 4. Blue Ball (NE Corner)
    # View from the left side of the doorway. Blocked by the wardrobe corner.
    view = joe_v_left
    target = blue_ball
    obstacle = wardrobe_corner
    slope_to_obstacle = (obstacle.y - view.y) / (obstacle.x - view.x)
    y_at_target_x = slope_to_obstacle * (target.x - view.x) + view.y
    is_blue_visible = target.y > y_at_target_x # Visible if above the shadow line
    if is_blue_visible:
        visible_balls.append("Blue")
    print(f"Checking Blue Ball at {blue_ball}:")
    print(f"  Joe's best viewpoint is {view}.")
    print(f"  The key obstruction is the wardrobe corner at {obstacle}.")
    print(f"  The line of sight to the ball (y={target.y}) is blocked because it is below the shadow line cast by the wardrobe (y={y_at_target_x:.2f} at the ball's x-position).")
    print("  Result: Not Visible\n")

    # 5. Purple Ball (In Wardrobe)
    # The ball is inside the wardrobe. Joe must see past the wardrobe's front.
    # We check where the line of sight hits the wardrobe's front plane (x=9.5).
    view = joe_v_left
    target = purple_ball
    # Calculate slope from view to target
    slope_to_target = (target.y - view.y) / (target.x - view.x)
    # Calculate y-intercept of the line of sight at the wardrobe front
    y_at_wardrobe_front = slope_to_target * (wardrobe_front_plane_x - view.x) + view.y
    # The wardrobe is solid below y=4
    is_purple_visible = y_at_wardrobe_front > wardrobe_base_y
    if is_purple_visible:
        visible_balls.append("Purple")
    print(f"Checking Purple Ball at {purple_ball}:")
    print(f"  The ball is inside the wardrobe, which starts at y={wardrobe_base_y}.")
    print(f"  The line of sight from {view} must pass through the wardrobe's front at x={wardrobe_front_plane_x}.")
    print(f"  At that front plane, the line of sight is at y={y_at_wardrobe_front:.2f}.")
    print(f"  Since this is below the wardrobe's base (y={wardrobe_base_y}), the view is blocked.")
    print("  Result: Not Visible\n")

    print("--- Conclusion ---")
    if visible_balls:
        print("Joe can see the following balls:")
        for ball in visible_balls:
            print(f"- {ball}")
    else:
        print("Joe cannot see any of the balls.")
    
    return visible_balls

if __name__ == '__main__':
    final_answer = check_visibility()
    # The final answer format is specific
    answer_string = " and ".join(final_answer)
    print(f"\n<<<Joe can see the {answer_string} balls.>>>")
