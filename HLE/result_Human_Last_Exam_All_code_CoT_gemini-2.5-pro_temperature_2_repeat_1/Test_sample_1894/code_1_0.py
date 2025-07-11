import math

def analyze_visibility():
    """
    Analyzes which balls Joe can see based on the room's layout and obstacles.
    All dimensions are in feet.
    """

    # --- Setup of the Room and Objects ---
    # Joe's eye position: He can move in the 3ft doorway (4.5 to 7.5 on x-axis)
    # and lean in 3 inches (0.25 ft). His height is 5ft, eye level ~4.7ft.
    joe_y = 0.25
    joe_x_range = (4.5, 7.5)
    joe_z = 4.7
    vantage_points = {
        "left": (joe_x_range[0], joe_y, joe_z),
        "right": (joe_x_range[1], joe_y, joe_z),
    }

    # Center of balls (x, y, z)
    balls = {
        "Red": {"pos": (11.75, 0.25, 0.25), "visible": False, "reason": "Not yet analyzed."},
        "Blue": {"pos": (11.75, 11.75, 2.25), "visible": False, "reason": "Not yet analyzed."},
        "Yellow": {"pos": (0.25, 11.75, 0.25), "visible": False, "reason": "Not yet analyzed."},
        "Green": {"pos": (0.25, 0.25, 7.25), "visible": False, "reason": "Not yet analyzed."},
        "Purple": {"pos": (11.75, 4.25, 0.25), "visible": False, "reason": "Not yet analyzed."},
    }

    print("Analyzing visibility for Joe from the doorway...\n")
    visible_balls = []

    # --- Analysis for each ball ---

    # 1. Green Ball
    ball_name = "Green"
    ball_z = balls[ball_name]["pos"][2]
    bookshelf_height = 7.0
    print(f"--- Checking the {ball_name} ball ---")
    print(f"The {ball_name} ball is on top of a bookshelf of height {bookshelf_height:.1f} feet, placing the ball's center at z={ball_z:.2f} feet.")
    print(f"Joe's eye level is at z={joe_z:.2f} feet.")
    balls[ball_name]["reason"] = f"The ball at z={ball_z:.2f} is higher than Joe's eyes (z={joe_z:.2f}) and is on top of the 7-foot shelf, so he cannot see its top surface from his position."
    print("Conclusion: " + balls[ball_name]["reason"] + "\n")

    # 2. Red Ball
    ball_name = "Red"
    ball_pos = balls[ball_name]["pos"]
    print(f"--- Checking the {ball_name} ball ---")
    print(f"The {ball_name} ball is at ({ball_pos[0]:.2f}, {ball_pos[1]:.2f}) on the floor in the southeast corner.")
    print(f"Joe's eye position is at y={joe_y:.2f} and the ball is at y={ball_pos[1]:.2f}.")
    print("This creates a direct, straight line of sight along the south wall.")
    balls[ball_name]["reason"] = "No obstacles are between Joe and the corner."
    balls[ball_name]["visible"] = True
    print("Conclusion: The Red ball is visible.\n")

    # 3. Yellow Ball
    ball_name = "Yellow"
    ball_pos = balls[ball_name]["pos"]
    eye_pos = vantage_points["right"] # Most advantageous position
    print(f"--- Checking the {ball_name} ball ---")
    print(f"The {ball_name} ball is in the northwest corner at ({ball_pos[0]:.2f}, {ball_pos[1]:.2f}).")
    print(f"To see around the main door, Joe's best vantage point is the far right of the doorway at ({eye_pos[0]:.2f}, {eye_pos[1]:.2f}).")
    # Check intersection with main door's plane (x=4.5)
    m = (ball_pos[1] - eye_pos[1]) / (ball_pos[0] - eye_pos[0])
    x_intersect = 4.5
    y_intersect = m * (x_intersect - eye_pos[0]) + eye_pos[1]
    print(f"The line of sight from ({eye_pos[0]:.2f}, {eye_pos[1]:.2f}) to ({ball_pos[0]:.2f}, {ball_pos[1]:.2f}) is calculated.")
    print(f"Equation: y - {eye_pos[1]:.2f} = {m:.3f} * (x - {eye_pos[0]:.2f})")
    print(f"At the main door's plane (x={x_intersect:.1f}), the sight line passes at y={y_intersect:.2f}.")
    print("Since the door obstacle only extends to y=3.0, the view is clear past the door.")
    # Check if view is blocked by bookshelf (7ft high)
    mz = (ball_pos[2] - eye_pos[2]) / math.sqrt((ball_pos[0]-eye_pos[0])**2 + (ball_pos[1]-eye_pos[1])**2)
    # Check height at front of bookshelf x=1
    dist_to_shelf_plane = eye_pos[0] - 1.0
    dist_along_ray = dist_to_shelf_plane / ( (eye_pos[0]-ball_pos[0])/math.sqrt((ball_pos[0]-eye_pos[0])**2 + (ball_pos[1]-eye_pos[1])**2) )
    z_at_shelf = eye_pos[2] + mz * dist_along_ray
    print(f"The bookshelf is 7ft tall, but the line of sight passes well over it (at a height of z={z_at_shelf:.2f} ft at the bookshelf's x=1.0 plane).")
    balls[ball_name]["reason"] = "A clear line of sight exists from the right side of the doorway, passing over the main door and bookshelf."
    balls[ball_name]["visible"] = True
    print("Conclusion: The Yellow ball is visible.\n")

    # 4. Blue Ball
    ball_name = "Blue"
    ball_pos = balls[ball_name]["pos"]
    eye_pos = vantage_points["left"] # Most advantageous position
    print(f"--- Checking the {ball_name} ball ---")
    print(f"The {ball_name} ball is in the northeast corner at ({ball_pos[0]:.2f}, {ball_pos[1]:.2f}).")
    print(f"Joe's best vantage point is the far left of the doorway at ({eye_pos[0]:.2f}, {eye_pos[1]:.2f}).")
    # Check intersection with wardrobe's front plane (x=9.5)
    m = (ball_pos[1] - eye_pos[1]) / (ball_pos[0] - eye_pos[0])
    x_intersect = 9.5
    y_intersect = m * (x_intersect - eye_pos[0]) + eye_pos[1]
    print(f"The line of sight from ({eye_pos[0]:.2f}, {eye_pos[1]:.2f}) to ({ball_pos[0]:.2f}, {ball_pos[1]:.2f}) is calculated.")
    print(f"Equation: y - {eye_pos[1]:.2f} = {m:.3f} * (x - {eye_pos[0]:.2f})")
    print(f"At the wardrobe's front (x={x_intersect:.1f}), the sight line passes at y={y_intersect:.2f}.")
    print("Since the wardrobe body only extends to y=8.0, the view is clear, just skimming past the top-north corner of the wardrobe.")
    balls[ball_name]["reason"] = "A clear line of sight exists from the left side of the doorway."
    balls[ball_name]["visible"] = True
    print("Conclusion: The Blue ball is visible.\n")
    
    # 5. Purple Ball
    ball_name = "Purple"
    ball_pos = balls[ball_name]["pos"]
    eye_pos = vantage_points["left"]
    print(f"--- Checking the {ball_name} ball ---")
    print(f"The {ball_name} ball is inside the wardrobe at ({ball_pos[0]:.2f}, {ball_pos[1]:.2f}).")
    print("To see it, Joe must look through the wardrobe's opening, which is at x=7.5.")
    print("The southern wardrobe door forms a wall from x=7.5 to x=9.5 at y=4.0.")
    # Check intersection with the wardrobe opening plane (x=7.5)
    m = (ball_pos[1] - eye_pos[1]) / (ball_pos[0] - eye_pos[0])
    x_intersect = 7.5
    y_intersect = m * (x_intersect - eye_pos[0]) + eye_pos[1]
    print(f"The line of sight from the leftmost vantage point ({eye_pos[0]:.2f}, {eye_pos[1]:.2f}) is calculated.")
    print(f"Equation: y - {eye_pos[1]:.2f} = {m:.3f} * (x - {eye_pos[0]:.2f})")
    print(f"At the wardrobe opening (x={x_intersect:.1f}), the line of sight passes at y={y_intersect:.2f}.")
    balls[ball_name]["reason"] = f"This sight line (at y={y_intersect:.2f}) is below the opening (which starts at y=4.0), so it is blocked by the wardrobe's open door."
    print("Conclusion: " + balls[ball_name]["reason"] + "\n")

    # --- Final Summary ---
    print("="*30)
    print("Final Result:")
    for name, data in balls.items():
        if data["visible"]:
            visible_balls.append(name)
    
    if visible_balls:
        print(f"Joe can see the following balls: {', '.join(visible_balls)}.")
    else:
        print("Joe cannot see any of the balls.")
    
    return visible_balls

final_answer = analyze_visibility()
# The required output format.
print(f"\n<<<The visible balls are {', '.join(final_answer)}>>>")
