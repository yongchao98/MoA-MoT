import math

def check_visibility():
    """
    Analyzes and prints the visibility of each ball based on the room layout.
    """
    # --- Define Room and Objects (all units in feet) ---

    # Joe's extreme viewpoints in the doorway (x, y)
    joe_viewpoints = {
        "left_doorway": (4.5, 0.25),
        "right_doorway": (6.0, 0.25)
    }
    joe_eye_height = 5.0

    # Balls' locations (x, y)
    balls = {
        "Red": (11.75, 0.25),
        "Blue": (11.75, 11.75),
        "Yellow": (0.25, 11.75),
        "Green": (0.75, 0.25), # Approx. location on the shelf
        "Purple": (11.75, 4.25)
    }

    # Obstacles' key features
    bookshelf = {"x_front": 1.0, "y_min": 0, "y_max": 4.0, "z_height": 7.0}
    wardrobe = {"x_front": 9.5, "y_min": 4.0, "y_max": 8.0}

    print("Analyzing visibility for each ball...\n")

    visible_balls = []

    # --- Green Ball Analysis (Height Check) ---
    ball_name = "Green"
    print(f"--- Checking the {ball_name} ball ---")
    green_ball_shelf_height = bookshelf["z_height"]
    print(f"The {ball_name} ball is on a shelf {green_ball_shelf_height} feet tall.")
    print(f"Joe's eye height is {joe_eye_height} feet.")
    print(f"Equation: Is Joe's eye height ({joe_eye_height}) > shelf height ({green_ball_shelf_height})?")
    if joe_eye_height > green_ball_shelf_height:
        print("Result: Yes. Visibility is possible.")
    else:
        print(f"Result: No. Since {joe_eye_height} < {green_ball_shelf_height}, Joe cannot see over the shelf to an object on top of it.")
    print("Conclusion: The Green ball is NOT visible.\n")


    # --- Red Ball Analysis ---
    ball_name = "Red"
    ball_pos = balls[ball_name]
    viewpoint = joe_viewpoints["right_doorway"]
    print(f"--- Checking the {ball_name} ball ---")
    print(f"Target: {ball_name} ball at {ball_pos}")
    print(f"Viewpoint: Right side of doorway at {viewpoint}")
    print("The line of sight is a straight line at a constant height.")
    print(f"Equation: y = {ball_pos[1]}")
    print(f"The nearest obstacle is the wardrobe, which starts at y = {wardrobe['y_min']}.")
    print(f"Result: The line of sight y-coordinate ({ball_pos[1]}) is less than the wardrobe's minimum y-coordinate ({wardrobe['y_min']}).")
    print(f"Conclusion: The {ball_name} ball is VISIBLE.\n")
    visible_balls.append(ball_name)


    # --- Yellow Ball Analysis ---
    ball_name = "Yellow"
    ball_pos = balls[ball_name]
    viewpoint = joe_viewpoints["right_doorway"]
    obstacle = bookshelf
    print(f"--- Checking the {ball_name} ball ---")
    print(f"Target: {ball_name} ball at {ball_pos}")
    print(f"Viewpoint: Right side of doorway at {viewpoint}")
    print(f"The main potential obstacle is the bookshelf, with its front at x = {obstacle['x_front']}.")

    # Equation of the line of sight
    vx, vy = viewpoint
    tx, ty = ball_pos
    m = (ty - vy) / (tx - vx)
    
    # Calculate height of LOS at the obstacle's x-position
    x_obs = obstacle['x_front']
    y_los = m * (x_obs - vx) + vy

    print(f"Line of Sight Slope (m) = ({ty} - {vy}) / ({tx} - {vx}) = {m:.4f}")
    print(f"Equation for line-of-sight height (y) at the bookshelf (x={x_obs}):")
    print(f"y = {m:.4f} * ({x_obs} - {vx}) + {vy} = {y_los:.4f}")
    print(f"Result: The line of sight height ({y_los:.4f}) is greater than the bookshelf's max height ({obstacle['y_max']}). The view is clear.")
    print(f"Conclusion: The {ball_name} ball is VISIBLE.\n")
    visible_balls.append(ball_name)


    # --- Blue Ball Analysis ---
    ball_name = "Blue"
    ball_pos = balls[ball_name]
    viewpoint = joe_viewpoints["left_doorway"]
    obstacle = wardrobe
    print(f"--- Checking the {ball_name} ball ---")
    print(f"Target: {ball_name} ball at {ball_pos}")
    print(f"Viewpoint: Left side of doorway at {viewpoint}")
    print(f"The main potential obstacle is the wardrobe, with its top edge at y = {obstacle['y_max']} at x = {obstacle['x_front']}.")
    
    vx, vy = viewpoint
    tx, ty = ball_pos
    m = (ty - vy) / (tx - vx)
    
    x_obs = obstacle['x_front']
    y_los = m * (x_obs - vx) + vy
    
    print(f"Line of Sight Slope (m) = ({ty} - {vy}) / ({tx} - {vx}) = {m:.4f}")
    print(f"Equation for line-of-sight height (y) at the wardrobe front (x={x_obs}):")
    print(f"y = {m:.4f} * ({x_obs} - {vx}) + {vy} = {y_los:.4f}")
    print(f"Result: The line of sight height ({y_los:.4f}) is greater than the wardrobe's top edge ({obstacle['y_max']}). The view passes over the wardrobe.")
    print(f"Conclusion: The {ball_name} ball is VISIBLE.\n")
    visible_balls.append(ball_name)


    # --- Purple Ball Analysis ---
    ball_name = "Purple"
    ball_pos = balls[ball_name]
    viewpoint = joe_viewpoints["right_doorway"]
    obstacle = wardrobe
    print(f"--- Checking the {ball_name} ball ---")
    print(f"Target: {ball_name} ball at {ball_pos}, inside the wardrobe.")
    print(f"Viewpoint: Right side of doorway at {viewpoint}")
    print(f"The view is blocked if the line of sight passes below the wardrobe opening, which starts at y = {obstacle['y_min']} at x = {obstacle['x_front']}.")

    vx, vy = viewpoint
    tx, ty = ball_pos
    m = (ty - vy) / (tx - vx)

    x_obs = obstacle['x_front']
    y_los = m * (x_obs - vx) + vy

    print(f"Line of Sight Slope (m) = ({ty} - {vy}) / ({tx} - {vx}) = {m:.4f}")
    print(f"Equation for line-of-sight height (y) at the wardrobe front (x={x_obs}):")
    print(f"y = {m:.4f} * ({x_obs} - {vx}) + {vy} = {y_los:.4f}")
    print(f"Result: The line of sight height ({y_los:.4f}) is less than the wardrobe opening's bottom edge ({obstacle['y_min']}). The view is blocked by the wardrobe's frame.")
    print(f"Conclusion: The {ball_name} ball is NOT visible.\n")

    # --- Final Summary ---
    print("--- Summary ---")
    print("Based on the analysis, Joe can see the following balls:")
    for ball in visible_balls:
        print(f"- {ball}")

if __name__ == '__main__':
    check_visibility()