import math

def check_visibility():
    """
    Analyzes the room layout to determine which balls are visible to Joe.
    Prints the step-by-step analysis and the final result.
    """
    # Define key locations using a coordinate system where (0,0) is the SW corner.
    # Joe's viewpoint is a line segment from x=4.5 to x=7.5 at y=-0.25 (3" into the room)
    joe_view_west = (4.5, -0.25)
    joe_view_east = (7.5, -0.25)
    joe_eye_height = 5.0

    # Obstacle definitions
    # Bookshelf is 7ft tall, taller than Joe's eyes.
    bookshelf_height = 7.0
    # Wardrobe doors create an opening on the line x=7.5, between y=4 and y=8
    wardrobe_opening_y = (4.0, 8.0)
    wardrobe_corner = (9.5, 8.0) # Top-left corner of the wardrobe body
    # Main door opens inward from x=4.5, creating a barrier up to y=3
    door_block = {'x': 4.5, 'y_max': 3.0}

    # Ball definitions (name, (x, y), height_z)
    balls = {
        "Green": ((0.25, 0.25), 7.25),
        "Purple": ((11.75, 4.25), 0.25),
        "Red": ((11.75, 0.25), 0.25),
        "Yellow": ((0.25, 11.75), 0.25),
        "Blue": ((11.75, 11.75), 2.75)
    }

    visible_balls = []

    print("Analyzing visibility for each ball:\n")

    # --- Green Ball ---
    green_ball_pos, green_ball_height = balls["Green"]
    print("1. Green Ball (on the 7ft bookshelf):")
    # Vertical check
    if joe_eye_height < bookshelf_height:
        print(f"   Joe's eye height ({joe_eye_height} ft) is below the top of the bookshelf ({bookshelf_height} ft).")
        print("   The bookshelf itself blocks the line of sight to the green ball placed on top of it.")
        print("   => Green Ball is NOT visible.\n")
    else: # This branch is not taken, but included for completeness
        visible_balls.append("Green")

    # --- Purple Ball ---
    purple_ball_pos, _ = balls["Purple"]
    print("2. Purple Ball (inside the wardrobe):")
    # To see the ball, Joe must see past the edge of the wardrobe door opening at (7.5, 4.0).
    # We check from Joe's most advantageous viewpoint for this, which is the far west side of the doorway.
    viewer = joe_view_west
    target = purple_ball_pos
    # Calculate line of sight
    m = (target[1] - viewer[1]) / (target[0] - viewer[0])
    # y = m*(x - x1) + y1. Find y at the wardrobe opening's x-coordinate.
    y_at_opening = m * (7.5 - viewer[0]) + viewer[1]
    print(f"   The line of sight from Joe's westernmost viewpoint ({viewer}) must pass through the wardrobe opening at x=7.5.")
    print(f"   At x=7.5, the line of sight has a y-coordinate of {y_at_opening:.2f}.")
    print(f"   The wardrobe opening begins at y={wardrobe_opening_y[0]:.2f}.")
    if y_at_opening < wardrobe_opening_y[0]:
        print("   Since the line of sight is below the opening, it is blocked by the wardrobe's lower door.")
        print("   => Purple Ball is NOT visible.\n")
    else: # Not taken
        visible_balls.append("Purple")
        
    # --- Red Ball ---
    print("3. Red Ball (in the southeast corner):")
    print("   The line of sight is from the doorway (~y=0) to the corner (~y=0).")
    print("   The main obstacles are the bookshelf (x<1), the open door (x=4.5), and the wardrobe (y>4).")
    print("   None of these obstacles intersect the low-angle line of sight to the southeast corner.")
    print("   => Red Ball is VISIBLE.\n")
    visible_balls.append("Red")

    # --- Yellow Ball ---
    yellow_ball_pos, _ = balls["Yellow"]
    print("4. Yellow Ball (in the northwest corner):")
    # To see this, Joe needs to look past the open door at x=4.5. His best view is from the east.
    viewer = joe_view_east
    target = yellow_ball_pos
    m = (target[1] - viewer[1]) / (target[0] - viewer[0])
    # y = m*(x - x1) + y1. Find y at the door's x-coordinate.
    y_at_door = m * (door_block['x'] - viewer[0]) + viewer[1]
    print(f"   From Joe's easternmost viewpoint ({viewer}), we check if the view is blocked by the open door at x={door_block['x']}.")
    print(f"   At x={door_block['x']}, the line of sight has a y-coordinate of {y_at_door:.2f}.")
    print(f"   The door only blocks up to y={door_block['y_max']:.2f}.")
    if y_at_door > door_block['y_max']:
        print("   The line of sight passes above the door. The bookshelf (x<1) is also not in the way.")
        print("   => Yellow Ball is VISIBLE.\n")
        visible_balls.append("Yellow")
    else: # Not taken
        print("   => Yellow Ball is NOT visible.\n")
        
    # --- Blue Ball ---
    blue_ball_pos, _ = balls["Blue"]
    print("5. Blue Ball (in the northeast corner):")
    # To see this, Joe must look through the wardrobe opening. His best view is from the west.
    viewer = joe_view_west
    target = blue_ball_pos
    
    # Check 1: Does it pass through the wardrobe opening (x=7.5, y=[4,8])?
    m_ball = (target[1] - viewer[1]) / (target[0] - viewer[0])
    y_at_opening = m_ball * (7.5 - viewer[0]) + viewer[1]
    print(f"   From Joe's westernmost viewpoint ({viewer}), we check if the view is blocked by the wardrobe.")
    print(f"   At the wardrobe opening (x=7.5), the line of sight has a y-coordinate of {y_at_opening:.2f}.")
    
    if wardrobe_opening_y[0] < y_at_opening < wardrobe_opening_y[1]:
        print(f"   This is within the opening's range of y={wardrobe_opening_y}. The view is not blocked by the doors.")
        # Check 2: Does it clear the wardrobe body corner at (9.5, 8)?
        slope_to_corner = (wardrobe_corner[1] - viewer[1]) / (wardrobe_corner[0] - viewer[0])
        print(f"   Now we check if the line of sight clears the wardrobe's top-left corner at {wardrobe_corner}.")
        print(f"   Slope to ball = {m_ball:.4f}. Slope to corner = {slope_to_corner:.4f}.")
        if m_ball > slope_to_corner:
            print("   The line to the ball is steeper than the line to the corner, so it passes above it.")
            print("   => Blue Ball is VISIBLE.\n")
            visible_balls.append("Blue")
        else: # Not taken
            print("   The line of sight is blocked by the wardrobe body.\n")
            print("   => Blue Ball is NOT visible.\n")
    else: # Not taken
        print(f"   This is outside the opening's range of y={wardrobe_opening_y}. The view is blocked by a door.\n")
        print("   => Blue Ball is NOT visible.\n")
        
    # --- Final Result ---
    print("--------------------------------------------------")
    print("Final Answer: The balls Joe can see are:")
    for ball in visible_balls:
        print(f"- {ball}")
    print("--------------------------------------------------")
    
    return visible_balls

if __name__ == '__main__':
    final_answer = check_visibility()
    # The final output format requested by the user
    # Convert list to a comma-separated string for the final answer block
    answer_string = ", ".join(sorted(final_answer))
    print(f"<<<{answer_string}>>>")
