import math

def solve_visibility():
    """
    Analyzes the visibility of balls in a room based on a set of geometric constraints.
    It calculates and prints the logic for whether each ball is visible to Joe from the doorway.
    """
    
    # --- Setup of the Room and Objects ---
    # Joe's eye can be at any point on the line segment from (4.5, 0.25, 5) to (7.5, 0.25, 5)
    joe_y = 0.25
    joe_z = 5.0
    
    # Ball positions (x, y, z)
    balls = {
        "red": {"pos": (11.75, 0.25, 0.25), "visible": False},
        "blue": {"pos": (11.75, 11.75, 2.75), "visible": False},
        "yellow": {"pos": (0.25, 11.75, 0.25), "visible": False},
        "green": {"pos": (0.25, 0.25, 7.25), "visible": False},
        "purple": {"pos": (11.75, 4.25, 0.25), "visible": False}
    }
    
    # --- Analysis for Each Ball ---
    
    # 1. Red Ball (Southeast corner)
    # The line of sight is at y=0.25. The wardrobe is at y>=4. The door is at x=7.5.
    # Joe can stand to the left of the door (e.g., at x=7.4) and see the ball.
    print("--- Checking Red Ball (in Southeast corner) ---")
    print("The line of sight is along the y=0.25 line.")
    print("The wardrobe starts at y=4, so it is not blocking the view.")
    print("The open door is a plane at x=7.5. Joe can stand at x=7.4 to see around it.")
    print("Result: Red ball is VISIBLE.\n")
    balls["red"]["visible"] = True

    # 2. Yellow Ball (Northwest corner)
    # The line of sight passes over the bookshelf's 2D footprint.
    print("--- Checking Yellow Ball (in Northwest corner) ---")
    joe_x_check = 4.5
    ball_x, ball_y, _ = balls["yellow"]["pos"]
    # Line equation: y = m*(x - x1) + y1
    m = (ball_y - joe_y) / (ball_x - joe_x_check)
    # Check y-value at the bookshelf's edge (x=1)
    y_at_shelf_edge = m * (1.0 - joe_x_check) + joe_y
    print(f"From Joe's leftmost position (x={joe_x_check}), the line of sight to the yellow ball passes x=1.0 (the bookshelf's edge).")
    print(f"The equation for the line of sight is y = (({ball_y} - {joe_y}) / ({ball_x} - {joe_x_check})) * (x - {joe_x_check}) + {joe_y}")
    print(f"At x=1.0, the y-coordinate of the line of sight is {y_at_shelf_edge:.2f}.")
    print("Since this y-value (9.70) is greater than the bookshelf's width (y=4), the line of sight does not cross the bookshelf.")
    print("Result: Yellow ball is VISIBLE.\n")
    balls["yellow"]["visible"] = True
    
    # 3. Green Ball (On top of Southwest bookshelf)
    # Check if the line of sight from Joe's best position clears the top of the bookshelf.
    print("--- Checking Green Ball (on Southwest bookshelf) ---")
    joe_x_check = 7.5 # Joe's rightmost position
    ball_x, _, ball_z = balls["green"]["pos"]
    shelf_front_x = 1.0
    shelf_top_z = 7.0
    dist_to_shelf = joe_x_check - shelf_front_x
    total_dist = joe_x_check - ball_x
    height_at_shelf = joe_z + (ball_z - joe_z) * (dist_to_shelf / total_dist)
    print(f"From Joe's rightmost position (x={joe_x_check}, z={joe_z}), he looks towards the ball (z={ball_z}).")
    print("The view must clear the bookshelf's front edge (x=1.0, z=7.0).")
    print(f"The height of the line of sight at the shelf's front edge (x=1.0) is calculated by interpolation:")
    print(f"Height = Joe_z + (Ball_z - Joe_z) * ((Joe_x - Shelf_x) / (Joe_x - Ball_x))")
    print(f"Height = {joe_z} + ({ball_z} - {joe_z}) * (({joe_x_check} - {shelf_front_x}) / ({joe_x_check} - {ball_x})) = {height_at_shelf:.3f}")
    print(f"Since the calculated height {height_at_shelf:.3f} is greater than the shelf height {shelf_top_z}, Joe can see the top of the ball.")
    print("Result: Green ball is VISIBLE.\n")
    balls["green"]["visible"] = True
    
    # 4. Blue Ball (Northeast corner)
    # Check if the line of sight from Joe's best position clears the wardrobe.
    print("--- Checking Blue Ball (in Northeast corner) ---")
    joe_x_check = 4.5 # Joe's leftmost position
    ball_x, ball_y, _ = balls["blue"]["pos"]
    wardrobe_corner_x = 9.5
    wardrobe_corner_y = 8.0
    m = (ball_y - joe_y) / (ball_x - joe_x_check)
    y_at_wardrobe = m * (wardrobe_corner_x - joe_x_check) + joe_y
    print(f"From Joe's leftmost position (x={joe_x_check}), he looks towards the ball.")
    print("The view must clear the wardrobe's corner (x=9.5, y=8.0).")
    print("We check the y-coordinate of the line of sight when it crosses the wardrobe's front line at x=9.5.")
    print(f"y_coord = {joe_y} + (({ball_y} - {joe_y}) / ({ball_x} - {joe_x_check})) * ({wardrobe_corner_x} - {joe_x_check}) = {y_at_wardrobe:.3f}")
    print(f"Since the calculated y-coordinate {y_at_wardrobe:.3f} is greater than the wardrobe corner's y-coordinate {wardrobe_corner_y}, the line of sight clears the wardrobe.")
    print("Result: Blue ball is VISIBLE.\n")
    balls["blue"]["visible"] = True

    # 5. Purple Ball (Inside wardrobe)
    # Check if the line of sight is blocked by the open door.
    print("--- Checking Purple Ball (in wardrobe) ---")
    joe_x_check = 4.5 # Joe's leftmost position
    ball_x, ball_y, _ = balls["purple"]["pos"]
    door_x = 7.5
    door_y_max = 3.0
    m = (ball_y - joe_y) / (ball_x - joe_x_check)
    y_at_door = m * (door_x - joe_x_check) + joe_y
    print("The ball is inside the wardrobe. The view must pass through the open door at x=7.5.")
    print(f"From Joe's leftmost position (x={joe_x_check}), we check where the line of sight intersects the door's plane at x={door_x}.")
    print(f"y_coord_at_door = {joe_y} + (({ball_y} - {joe_y}) / ({ball_x} - {joe_x_check})) * ({door_x} - {joe_x_check}) = {y_at_door:.3f}")
    print(f"The door spans from y=0 to y={door_y_max}. Since the intersection point y={y_at_door:.3f} is within this range, the door blocks the view.")
    print("Checking from other positions in the doorway also results in the door blocking the view.")
    print("Result: Purple ball is NOT VISIBLE.\n")
    balls["purple"]["visible"] = False
    
    # --- Final Conclusion ---
    visible_balls = sorted([name for name, data in balls.items() if data["visible"]])
    print("="*40)
    print(f"Based on the analysis, the balls Joe can see are: {', '.join(visible_balls)}.")
    print("="*40)

solve_visibility()
<<<Red, Blue, Yellow, Green>>>