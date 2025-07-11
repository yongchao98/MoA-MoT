import math

def solve_visibility_problem():
    """
    Analyzes the visibility of balls in a room based on a set of rules.
    """
    # --- Setup: Define room, objects, and viewpoints ---

    # Joe's eye level and viewing position
    # Doorway is 3ft wide, centered on 12ft south wall (x=0 to 12).
    # Center is x=6. Doorway is from x=4.5 to x=7.5.
    # Joe can lean in 3 inches (0.25 ft).
    JOE_EYE_HEIGHT = 4.75  # in feet
    VIEWPOINT_WEST = (4.5, 0.25)  # Westernmost point Joe can see from
    VIEWPOINT_EAST = (7.5, 0.25)  # Easternmost point Joe can see from

    # Obstruction definitions (as points, lines, or polygons)
    # Bookshelf is 1ft deep, 4ft wide, 7ft tall, in SW corner.
    BOOKSHELF_HEIGHT = 7.0
    BOOKSHELF_CORNER = (1.0, 4.0)  # The corner that obstructs the view to the NW

    # Door is 3ft wide, hinged on east side of doorway (x=7.5), open 90 degrees.
    # It creates a 3ft visual barrier into the room.
    DOOR_EDGE = (7.5, 3.0) # The edge of the door obstructing view to the east

    # Wardrobe doors are open. The south door is the most likely obstruction for the NE view.
    # It's hinged at (9.5, 4) and opens 90 degrees, so its edge is at (7.5, 4).
    WARDROBE_DOOR_CORNER = (7.5, 4.0)
    
    # Ball locations
    BALLS = {
        "Green": {"pos": (0.25, 0.25), "location_desc": "on top of the 7-foot-tall bookshelf"},
        "Purple": {"pos": (11.75, 4.25), "location_desc": "in the southeast corner of the wardrobe"},
        "Red": {"pos": (11.75, 0.25), "location_desc": "in the southeast corner"},
        "Yellow": {"pos": (0.25, 11.75), "location_desc": "in the northwest corner"},
        "Blue": {"pos": (11.75, 11.75), "location_desc": "in the northeast corner"}
    }

    visible_balls = []

    print("Analyzing which balls Joe can see from the doorway.\n")

    # --- Analysis for each ball ---

    # 1. Green Ball
    print("--- Checking the Green ball ---")
    print(f"The Green ball is {BALLS['Green']['location_desc']}.")
    print(f"The bookshelf height is {BOOKSHELF_HEIGHT} feet.")
    print(f"Joe's eye level is {JOE_EYE_HEIGHT} feet.")
    if BALLS['Green']['location_desc'].startswith("on top") and BOOKSHELF_HEIGHT > JOE_EYE_HEIGHT:
        print(f"Since the bookshelf ({BOOKSHELF_HEIGHT} ft) is taller than Joe's eyes ({JOE_EYE_HEIGHT} ft), he cannot see the Green ball on top of it.")
        print("Conclusion: Green ball is NOT visible.\n")
    else:
        # This part would run if the ball was not on top, but we know it is.
        pass

    # 2. Purple Ball
    print("--- Checking the Purple ball ---")
    print(f"The Purple ball is {BALLS['Purple']['location_desc']}.")
    print("The wardrobe is an opaque object. Joe cannot see inside it from the doorway.")
    print("Conclusion: Purple ball is NOT visible.\n")

    # 3. Yellow Ball (in NW corner)
    print("--- Checking the Yellow ball ---")
    print(f"The Yellow ball is {BALLS['Yellow']['location_desc']}.")
    print("To see it, Joe must look past the bookshelf. His most advantageous viewpoint is the easternmost one.")
    V = VIEWPOINT_EAST
    T = BALLS['Yellow']['pos']
    P = BOOKSHELF_CORNER
    print(f"Viewpoint (V): {V}")
    print(f"Target (T): {T}")
    print(f"Obstruction Corner (P): {P}")
    print("Calculating the line of sight from V past P to the back wall (y=12).")
    
    # Equation of line from V to P: y - V_y = m * (x - V_x)
    m = (P[1] - V[1]) / (P[0] - V[0])
    print(f"The slope (m) of the line from V to P is ({P[1]} - {V[1]}) / ({P[0]} - {V[0]}) = {m:.3f}")
    
    # Find x where this line intersects y=12 (the back wall)
    # 12 - V_y = m * (x - V_x) => x = (12 - V_y)/m + V_x
    x_intercept = (12 - V[1]) / m + V[0]
    print(f"The equation for the line of sight grazing the obstacle is: y - {V[1]} = {m:.3f} * (x - {V[0]})")
    print(f"This sight line intersects the back wall (y=12) at x = {x_intercept:.2f}.")
    print(f"The Yellow ball is at x={T[0]}.")
    if T[0] > x_intercept:
        print(f"Since {T[0]} > {x_intercept:.2f}, the ball is not in the bookshelf's shadow.")
        print("Conclusion: Yellow ball is VISIBLE.\n")
        visible_balls.append("Yellow")
    else:
        print(f"Since {T[0]} <= {x_intercept:.2f}, the ball is in the bookshelf's shadow.")
        print("Conclusion: Yellow ball is NOT visible.\n")
        
    # 4. Red Ball (in SE corner)
    print("--- Checking the Red ball ---")
    print(f"The Red ball is {BALLS['Red']['location_desc']}.")
    print("To see it, Joe must look past the open door. His most advantageous viewpoint is the westernmost one.")
    V = VIEWPOINT_WEST
    T = BALLS['Red']['pos']
    P = DOOR_EDGE
    print(f"Viewpoint (V): {V}")
    print(f"Target (T): {T}")
    print(f"Obstruction Corner (P): {P}")
    print("Calculating the line of sight from V past P.")
    
    # Equation of line from V to P: y - V_y = m * (x - V_x)
    m = (P[1] - V[1]) / (P[0] - V[0])
    print(f"The slope (m) of the line from V to P is ({P[1]} - {V[1]}) / ({P[0]} - {V[0]}) = {m:.3f}")
    print(f"The equation for the line of sight grazing the obstacle is: y - {V[1]} = {m:.3f} * (x - {V[0]})")

    # Find the y-value of this sight line at the ball's x-coordinate
    y_at_target_x = m * (T[0] - V[0]) + V[1]
    print(f"At the ball's x-coordinate ({T[0]}), this sight line is at y = {y_at_target_x:.2f}.")
    print(f"The Red ball is at y={T[1]}.")
    if T[1] < y_at_target_x:
        print(f"Since the ball's y-coordinate ({T[1]}) is less than the sight line's y-coordinate ({y_at_target_x:.2f}), the ball is blocked by the door.")
        print("Conclusion: Red ball is NOT visible.\n")
    else:
        print("Conclusion: Red ball is VISIBLE.\n")
        visible_balls.append("Red")
        
    # 5. Blue Ball (in NE corner)
    print("--- Checking the Blue ball ---")
    print(f"The Blue ball is {BALLS['Blue']['location_desc']}.")
    print("To see it, Joe must look past the open door and the wardrobe. His most advantageous viewpoint is the westernmost one.")
    V = VIEWPOINT_WEST
    T = BALLS['Blue']['pos']
    P = WARDROBE_DOOR_CORNER # The corner of the wardrobe's open door is the primary obstruction
    print(f"Viewpoint (V): {V}")
    print(f"Target (T): {T}")
    print(f"Obstruction Corner (P): {P}")
    print("We must check if the line of sight to the ball (T) is blocked by the obstruction (P).")
    
    # Calculate slope to target vs slope to obstruction corner
    slope_to_target = (T[1] - V[1]) / (T[0] - V[0])
    slope_to_obstruction = (P[1] - V[1]) / (P[0] - V[0])
    
    print(f"The slope of the line to the obstruction is m1 = ({P[1]} - {V[1]}) / ({P[0]} - {V[0]}) = {slope_to_obstruction:.3f}")
    print(f"The slope of the line to the target ball is m2 = ({T[1]} - {V[1]}) / ({T[0]} - {V[0]}) = {slope_to_target:.3f}")

    if slope_to_target > slope_to_obstruction:
        print(f"Since the slope to the ball ({slope_to_target:.3f}) is greater than the slope to the obstruction ({slope_to_obstruction:.3f}), the line of sight passes clear of the wardrobe.")
        print("Conclusion: Blue ball is VISIBLE.\n")
        visible_balls.append("Blue")
    else:
        print("Since the slope to the ball is not greater than the slope to the obstruction, the view is blocked.")
        print("Conclusion: Blue ball is NOT visible.\n")

    # --- Final Summary ---
    print("-------------------------------------")
    print("Final Result:")
    if visible_balls:
        print("Joe can see the following balls: " + ", ".join(visible_balls) + ".")
    else:
        print("Joe cannot see any of the balls.")
    print("-------------------------------------")

if __name__ == '__main__':
    solve_visibility_problem()