def solve_visibility():
    """
    Analyzes the visibility of several balls in a room from Joe's perspective in a doorway.
    This function models the room, objects, and Joe's position to calculate lines of sight.
    """
    # --- 1. Define Scene ---
    # Joe's eye level and viewing range in feet
    JOE_EYE_LEVEL = 4.8  # 5ft tall person, eyes slightly lower
    JOE_LEAN_IN = 0.25  # 3 inches
    JOE_VIEW_LEFT_2D = (4.5, JOE_LEAN_IN)  # x, y
    JOE_VIEW_RIGHT_2D = (7.5, JOE_LEAN_IN)  # x, y

    # Ball positions (x, y, z)
    # 6-inch diameter = 0.5 ft -> radius = 0.25 ft
    BALL_RADIUS = 0.25
    BALLS = {
        "red": {"pos": (12 - BALL_RADIUS, BALL_RADIUS, BALL_RADIUS), "desc": "in the southeast corner on the floor"},
        "blue": {"pos": (12 - BALL_RADIUS, 12 - BALL_RADIUS, 2.5 + BALL_RADIUS), "desc": "in the far northeast corner on a table"},
        "yellow": {"pos": (BALL_RADIUS, 12 - BALL_RADIUS, BALL_RADIUS), "desc": "in the northwest corner on the floor"},
        "green": {"pos": (BALL_RADIUS, BALL_RADIUS, 7 + BALL_RADIUS), "desc": "in the southwest corner on top of the 7ft shelf"},
        "purple": {"pos": (12 - BALL_RADIUS, 4 + BALL_RADIUS, BALL_RADIUS), "desc": "in the southeast corner of the wardrobe"}
    }
    
    # Obstacle key dimensions
    BOOKSHELF_FRONT_X = 1
    BOOKSHELF_TOP_Y = 4
    BOOKSHELF_TOP_Z = 7

    WARDROBE_DOOR_SOUTH_Y = 4
    WARDROBE_DOOR_SOUTH_X_MIN = 9.5
    WARDROBE_DOOR_SOUTH_X_MAX = 11.5

    WARDROBE_NORTH_CORNER_X = 9.5
    WARDROBE_NORTH_CORNER_Y = 8
    
    visible_balls = []

    print("Analyzing visibility of each ball from the doorway...")

    # --- Red Ball ---
    print("\n--- Checking Red ball ---")
    print(f"The red ball is {BALLS['red']['desc']} at y={BALLS['red']['pos'][1]:.2f}.")
    print(f"Joe can look straight ahead from the doorway (y={JOE_LEAN_IN}) to the ball.")
    print("There are no obstacles along this direct line of sight.")
    print("Result: Red ball is VISIBLE.")
    visible_balls.append("Red")

    # --- Green Ball ---
    print("\n--- Checking Green ball ---")
    print(f"The green ball is {BALLS['green']['desc']}.")
    print(f"The bookshelf height is {BOOKSHELF_TOP_Z} feet. Joe's eye level is {JOE_EYE_LEVEL} feet.")
    print(f"Since Joe's eyes ({JOE_EYE_LEVEL} ft) are below the top of the bookshelf ({BOOKSHELF_TOP_Z} ft), he cannot see an object placed on top of it.")
    print("Result: Green ball is NOT VISIBLE.")

    # --- Yellow Ball ---
    print("\n--- Checking Yellow ball ---")
    print(f"The yellow ball is {BALLS['yellow']['desc']}.")
    print("The main potential obstruction is the bookshelf in the southwest corner.")
    print("Checking line of sight from Joe's optimal (rightmost) position to the ball...")
    
    joe_pos = JOE_VIEW_RIGHT_2D
    ball_pos = (BALLS['yellow']['pos'][0], BALLS['yellow']['pos'][1])
    
    # Calculate y-value of the line of sight at the bookshelf's front edge (x=1)
    x1, y1 = joe_pos
    x2, y2 = ball_pos
    x = BOOKSHELF_FRONT_X
    m = (y2 - y1) / (x2 - x1)
    y_at_blocker = m * (x - x1) + y1
    
    print(f"Line of sight is from Joe at ({x1}, {y1}) to Ball at ({x2:.2f}, {y2:.2f}).")
    print(f"Equation: y - {y1} = (({y2:.2f} - {y1}) / ({x2:.2f} - {x1})) * (x - {x1})")
    print(f"At the bookshelf's front edge (x = {x}), we solve for y:")
    print(f"y = (({y2-y1:.2f}) / ({x2-x1:.2f})) * ({x} - {x1}) + {y1} = {y_at_blocker:.2f}")

    if y_at_blocker > BOOKSHELF_TOP_Y:
        print(f"Since the line of sight passes at y={y_at_blocker:.2f}, which is north of the bookshelf's edge at y={BOOKSHELF_TOP_Y}, the view is clear.")
        print("Result: Yellow ball is VISIBLE.")
        visible_balls.append("Yellow")
    else:
        print("Result: Yellow ball is NOT VISIBLE.")

    # --- Purple Ball ---
    print("\n--- Checking Purple ball ---")
    print(f"The purple ball is {BALLS['purple']['desc']}.")
    print("The main obstruction is the open southern wardrobe door, a vertical plane at y=4.")
    print("Checking line of sight from Joe's optimal (rightmost) position to the ball...")
    
    joe_pos = JOE_VIEW_RIGHT_2D
    ball_pos = (BALLS['purple']['pos'][0], BALLS['purple']['pos'][1])

    # Calculate x-value of the line of sight at the door's plane (y=4)
    x1, y1 = joe_pos
    x2, y2 = ball_pos
    y = WARDROBE_DOOR_SOUTH_Y
    m_inv = (x2 - x1) / (y2 - y1)
    x_at_blocker = m_inv * (y - y1) + x1

    print(f"Line of sight is from Joe at ({x1}, {y1}) to Ball at ({x2:.2f}, {y2:.2f}).")
    print(f"At the door's plane (y = {y}), we solve for x:")
    print(f"x = (({x2-x1:.2f}) / ({y2-y1:.2f})) * ({y} - {y1}) + {x1} = {x_at_blocker:.2f}")
    
    if WARDROBE_DOOR_SOUTH_X_MIN <= x_at_blocker <= WARDROBE_DOOR_SOUTH_X_MAX:
        print(f"The intersection at x={x_at_blocker:.2f} is within the door's span of [{WARDROBE_DOOR_SOUTH_X_MIN}, {WARDROBE_DOOR_SOUTH_X_MAX}].")
        print("Result: Purple ball is NOT VISIBLE.")
    else:
        print("Result: Purple ball is VISIBLE.") # This case should not be reached

    # --- Blue Ball ---
    print("\n--- Checking Blue ball ---")
    print(f"The blue ball is {BALLS['blue']['desc']}.")
    print("The main obstruction is the corner of the wardrobe and its open northern door.")
    print("Checking line of sight from Joe's optimal (leftmost) position to the ball...")

    joe_pos = JOE_VIEW_LEFT_2D
    ball_pos = (BALLS['blue']['pos'][0], BALLS['blue']['pos'][1])

    # Calculate y-value of the line of sight at the wardrobe's corner (x=9.5)
    x1, y1 = joe_pos
    x2, y2 = ball_pos
    x = WARDROBE_NORTH_CORNER_X
    m = (y2 - y1) / (x2 - x1)
    y_at_blocker = m * (x - x1) + y1

    print(f"Line of sight is from Joe at ({x1}, {y1}) to Ball at ({x2:.2f}, {y2:.2f}).")
    print(f"At the wardrobe's front corner (x = {x}), we solve for y:")
    print(f"y = (({y2-y1:.2f}) / ({x2-x1:.2f})) * ({x} - {x1}) + {y1} = {y_at_blocker:.2f}")
    
    if y_at_blocker > WARDROBE_NORTH_CORNER_Y:
        print(f"Since the line of sight at y={y_at_blocker:.2f} passes north of the wardrobe's corner at y={WARDROBE_NORTH_CORNER_Y}, the view is clear.")
        print("Result: Blue ball is VISIBLE.")
        visible_balls.append("Blue")
    else:
        print("Result: Blue ball is NOT VISIBLE.")

    # --- Final Result ---
    print("\n---------------------")
    print("SUMMARY")
    print("---------------------")
    if not visible_balls:
        print("Joe cannot see any of the balls.")
    else:
        final_answer = ", ".join(sorted(visible_balls))
        print(f"The balls Joe can see are: {final_answer}.")
        print(f"<<<{final_answer}>>>")

solve_visibility()