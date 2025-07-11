import math

def solve_visibility_puzzle():
    """
    Determines which balls Joe can see based on the room's layout and obstructions.
    The logic is explained step-by-step with calculations for each ball.
    """
    print("Analyzing the room to determine which balls Joe can see...\n")

    # --- Define Scene Coordinates ---
    # Joe's eye level is ~4.75 ft. He can move in the doorway [x=4.5 to 7.5]
    # and lean in 3 inches (y=0.25 ft).
    # We test from the most extreme valid viewpoints.
    joe_viewpoints = {
        "left": (4.5, 0.25),
        "right": (7.5, 0.25)
    }
    
    # Ball positions (x, y, z)
    balls = {
        "Red": {"pos": (11.75, 0.25, 0.25), "desc": "in the southeast corner on the floor"},
        "Blue": {"pos": (11.75, 11.75, 2.75), "desc": "in the northeast corner on a table"},
        "Yellow": {"pos": (0.25, 11.75, 0.25), "desc": "in the northwest corner on the floor"},
        "Green": {"pos": (0.25, 0.25, 7.25), "desc": "in the southwest corner on top of the shelf"},
        "Purple": {"pos": (11.75, 4.25, 0.25), "desc": "in the southeast corner of the wardrobe"}
    }

    # Obstruction definitions (2D top-down view)
    # Door is hinged at x=4.5 and opens inward, creating a wall from y=0 to y=3
    door_obstruction = {"x": 4.5, "y_range": (0, 3)}
    # Wardrobe door edge that may block the purple ball
    wardrobe_door_edge = (11.5, 4.0)
    # Wardrobe corner and door that may block the blue ball
    wardrobe_corner = (9.5, 8.0)


    visible_balls = []

    # --- Analysis for each ball ---

    # 1. Green Ball
    print("--- Checking Green Ball ---")
    print(f"Location: {balls['Green']['desc']} at {balls['Green']['pos']}")
    green_ball_pos = balls['Green']['pos']
    # To see the green ball at x=0.25, Joe at x>=4.5 must look through the door's plane at x=4.5.
    # We check from Joe's right-most viewpoint (7.5, 0.25) to the ball (0.25, 0.25).
    # The line of sight is a horizontal line at y=0.25.
    print("Logic: The line of sight is a horizontal line at y=0.25.")
    print(f"This line intersects the door's plane (at x={door_obstruction['x']}) at y=0.25.")
    print(f"The door obstruction exists from y={door_obstruction['y_range'][0]} to y={door_obstruction['y_range'][1]}.")
    if door_obstruction['y_range'][0] <= 0.25 <= door_obstruction['y_range'][1]:
        print("Result: The intersection point is within the door's bounds. The view is BLOCKED.")
    else:
        print("Result: The intersection point is outside the door's bounds. The view is CLEAR.")
        # This else case won't be hit based on the logic.
    print("-" * 25 + "\n")


    # 2. Yellow Ball
    print("--- Checking Yellow Ball ---")
    print(f"Location: {balls['Yellow']['desc']} at {balls['Yellow']['pos']}")
    yellow_ball_pos = balls['Yellow']['pos']
    joe_vp = joe_viewpoints["right"]
    # To see around the door, Joe should be at his right-most position.
    print(f"Logic: Check line of sight from Joe's right viewpoint {joe_vp} to the ball {yellow_ball_pos[:2]}.")
    print(f"The main obstruction is the door at x={door_obstruction['x']}.")
    print("Calculating the y-coordinate of the line of sight where it crosses the door's plane (x=4.5):")
    # Equation: y = m*(x - x1) + y1
    x1, y1 = joe_vp
    x2, y2 = yellow_ball_pos[:2]
    target_x = door_obstruction['x']
    y_intersect = ((y2 - y1) / (x2 - x1)) * (target_x - x1) + y1
    print(f"  y = (({y2} - {y1}) / ({x2} - {x1})) * ({target_x} - {x1}) + {y1}")
    print(f"Result: The line of sight crosses the door's plane at y = {y_intersect:.2f}.")
    if y_intersect > door_obstruction['y_range'][1]:
        print(f"This is higher than the door's edge at y={door_obstruction['y_range'][1]}. The view is CLEAR.")
        visible_balls.append("Yellow")
    else:
        print(f"This is within the door's range. The view is BLOCKED.")
    print("-" * 25 + "\n")


    # 3. Red Ball
    print("--- Checking Red Ball ---")
    print(f"Location: {balls['Red']['desc']} at {balls['Red']['pos']}")
    # The line of sight is along the south wall at y=0.25.
    # The wardrobe is at y=[4,8], so it does not block the view.
    print("Logic: The line of sight is along the y=0.25 line.")
    print("The nearest obstruction is the wardrobe, which starts at y=4.")
    print("Result: There are no obstructions in the path. The view is CLEAR.")
    visible_balls.append("Red")
    print("-" * 25 + "\n")


    # 4. Blue Ball
    print("--- Checking Blue Ball ---")
    print(f"Location: {balls['Blue']['desc']} at {balls['Blue']['pos']}")
    blue_ball_pos = balls['Blue']['pos']
    joe_vp = joe_viewpoints["left"]
    # To see around the wardrobe, Joe should be at his left-most position.
    print(f"Logic: Check line of sight from Joe's left viewpoint {joe_vp} to the ball {blue_ball_pos[:2]}.")
    print(f"The main obstruction is the wardrobe corner at {wardrobe_corner}.")
    print("Calculating the y-coordinate of the line of sight where it crosses the wardrobe's front (x=9.5):")
    x1, y1 = joe_vp
    x2, y2 = blue_ball_pos[:2]
    target_x = wardrobe_corner[0]
    y_intersect = ((y2 - y1) / (x2 - x1)) * (target_x - x1) + y1
    print(f"  y = (({y2} - {y1}) / ({x2} - {x1})) * ({target_x} - {x1}) + {y1}")
    print(f"Result: The line of sight crosses the wardrobe's front plane at y = {y_intersect:.2f}.")
    if y_intersect > wardrobe_corner[1]:
        print(f"This is higher than the wardrobe corner at y={wardrobe_corner[1]}. The view is CLEAR.")
        visible_balls.append("Blue")
    else:
        print(f"This is blocked by the wardrobe. The view is BLOCKED.")
    print("-" * 25 + "\n")


    # 5. Purple Ball
    print("--- Checking Purple Ball ---")
    print(f"Location: {balls['Purple']['desc']} at {balls['Purple']['pos']}")
    purple_ball_pos = balls['Purple']['pos']
    joe_vp = joe_viewpoints["right"]
    # To see into the wardrobe, Joe should be at his right-most position.
    print(f"Logic: Check line of sight from Joe's right viewpoint {joe_vp} to the ball {purple_ball_pos[:2]}.")
    print(f"The main obstruction is the open wardrobe door edge at {wardrobe_door_edge}.")
    print("Calculating the y-coordinate of the line of sight where it crosses the door edge's x-plane (x=11.5):")
    x1, y1 = joe_vp
    x2, y2 = purple_ball_pos[:2]
    target_x = wardrobe_door_edge[0]
    y_intersect = ((y2 - y1) / (x2 - x1)) * (target_x - x1) + y1
    print(f"  y = (({y2} - {y1}) / ({x2} - {x1})) * ({target_x} - {x1}) + {y1}")
    print(f"Result: The line of sight crosses the door's plane at y = {y_intersect:.2f}.")
    if y_intersect > wardrobe_door_edge[1]:
        print(f"This is higher than the door edge at y={wardrobe_door_edge[1]}. The view is CLEAR.")
        visible_balls.append("Purple")
    else:
        print(f"This is blocked by the door edge. The view is BLOCKED.")
    print("-" * 25 + "\n")


    # --- Final Conclusion ---
    print("="*40)
    print("Final Conclusion:")
    if visible_balls:
        print("Joe can see the following balls:")
        for ball in sorted(visible_balls):
            print(f"- {ball}")
    else:
        print("Joe cannot see any of the balls.")
    print("="*40)


if __name__ == '__main__':
    solve_visibility_puzzle()
    # The final answer is derived from the code's output.
    # Visible balls: Yellow, Red, Blue, Purple.
    final_answer = "Red, Blue, Yellow, and Purple"
    # print(f"\n<<<ANSWER>>>\n{final_answer}") # This is for the final output format.

<<<Joe can see the Red, Blue, Yellow, and Purple balls.>>>