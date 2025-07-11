import collections

# Define data structures for points and objects
Point = collections.namedtuple('Point', ['x', 'y', 'z'])
Ball = collections.namedtuple('Ball', ['name', 'color', 'pos'])
ViewRegion = collections.namedtuple('ViewRegion', ['x_min', 'x_max', 'y', 'z_eye'])

def analyze_visibility():
    """
    Models the room and calculates visibility of each ball from the doorway.
    """
    # --- Setup of the environment ---

    # Joe's viewpoint: in doorway [4.5, 7.5] on south wall (y=0), leans in 3" (0.25ft)
    # Eye height is 5ft person -> ~4.5ft
    joe_view = ViewRegion(x_min=4.5, x_max=7.5, y=0.25, z_eye=4.5)

    # Balls locations (center of 0.5ft diameter ball)
    balls = [
        Ball('Red', 'red', Point(11.75, 0.25, 0.25)),    # SE corner, on floor
        Ball('Green', 'green', Point(0.25, 0.25, 7.25)),  # SW corner, on 7ft shelf
        Ball('Yellow', 'yellow', Point(0.25, 11.75, 0.25)), # NW corner, on floor
        Ball('Blue', 'blue', Point(11.75, 11.75, 2.75)),  # NE corner, on table
        Ball('Purple', 'purple', Point(11.75, 4.25, 0.25)) # In SE corner of wardrobe
    ]
    
    visible_balls = []

    print("--- Analysis of Ball Visibility from the Doorway ---\n")
    print("ASSUMPTIONS:")
    print("- Room is 12x12ft. SW corner is origin (0,0).")
    print("- Joe's eye position (Viewpoint) can be at y=0.25ft, with x between 4.5ft and 7.5ft, at a height z=4.5ft.")
    print("- Door is 3ft wide, hinged at x=4.5, opens INWARD, forming a wall at x=4.5 from y=0 to y=3, height 6.7ft.")
    print("- Wardrobe (x=[9.5,12], y=[4,8]) has two doors hinged at the sides (y=4, y=8) that open 90 degrees outward, creating an opening at x=9.5.\n")

    # --- Ball by Ball Analysis ---

    # 1. Red Ball
    ball = balls[0]
    print(f"1. Checking {ball.color} Ball at ({ball.pos.x}, {ball.pos.y}, {ball.pos.z})")
    print("   The Red ball and Joe's eyes are at the same y-coordinate (y=0.25).")
    print("   The line of sight is a straight line along y=0.25.")
    print("   There are no obstructions (like the wardrobe) at this y-coordinate.")
    print(f"   Conclusion: {ball.color} ball is VISIBLE.\n")
    visible_balls.append(ball.color)

    # 2. Green Ball
    ball = balls[1]
    # View is blocked by the door at x=4.5. Best chance is to be at far right of doorway.
    view_pos = Point(joe_view.x_max, joe_view.y, joe_view.z_eye)
    obstruction_door_edge = Point(4.5, -1, 6.7) # x=4.5, z=6.7
    
    print(f"2. Checking {ball.color} Ball at ({ball.pos.x}, {ball.pos.y}, {ball.pos.z})")
    print(f"   To see to the left, Joe must look past the open door at x={obstruction_door_edge.x} ft.")
    print(f"   Best viewpoint is far right: ({view_pos.x}, {view_pos.y}, {view_pos.z}).")
    print("   We check if the line of sight passes over the door's top edge (z=6.7 ft).")
    
    vx, vz = view_pos.x, view_pos.z
    gx, gz = ball.pos.x, ball.pos.z
    ox = obstruction_door_edge.x
    oz = obstruction_door_edge.z

    # Line of sight equation in x-z plane: z = vz + (x - vx) * (gz - vz) / (gx - vx)
    sight_height_at_door = vz + (ox - vx) * (gz - vz) / (gx - vx)
    
    print(f"   Equation: z = {vz} + ({ox} - {vx}) * ({gz} - {vz}) / ({gx} - {vx})")
    print(f"   Calculated height of sightline at door (x={ox}): z = {sight_height_at_door:.2f} ft.")
    print(f"   The sightline height ({sight_height_at_door:.2f} ft) is lower than the door height ({oz} ft).")
    print(f"   Conclusion: {ball.color} ball is NOT VISIBLE (blocked by the door).\n")

    # 3. Yellow Ball
    ball = balls[2]
    # View is potentially blocked by bookshelf corner (x=1, y=4). Best viewpoint is far right.
    view_pos = Point(joe_view.x_max, joe_view.y, joe_view.z_eye)
    obstruction_bookshelf_corner = Point(1, 4, 7) # x=1, y=4
    
    print(f"3. Checking {ball.color} Ball at ({ball.pos.x}, {ball.pos.y}, {ball.pos.z})")
    print(f"   To see to the NW corner, the view might be blocked by the bookshelf corner at ({obstruction_bookshelf_corner.x}, {obstruction_bookshelf_corner.y}).")
    print(f"   Best viewpoint is far right: ({view_pos.x}, {view_pos.y}, {view_pos.z}).")
    
    vx, vy = view_pos.x, view_pos.y
    yx, yy = ball.pos.x, ball.pos.y
    ox, oy = obstruction_bookshelf_corner.x, obstruction_bookshelf_corner.y

    # Line of sight equation in x-y plane: y = vy + (x - vx) * (yy - vy) / (yx - vx)
    sight_y_at_obstacle = vy + (ox - vx) * (yy - vy) / (yx - vx)
    
    print(f"   Equation: y = {vy} + ({ox} - {vx}) * ({yy} - {vy}) / ({yx} - {vx})")
    print(f"   Calculated y-position of sightline at bookshelf's x-edge (x={ox}): y = {sight_y_at_obstacle:.2f} ft.")
    print(f"   The sightline y-position ({sight_y_at_obstacle:.2f} ft) is greater than the bookshelf's y-extent ({oy} ft).")
    print(f"   Conclusion: {ball.color} ball is VISIBLE.\n")
    visible_balls.append(ball.color)
    
    # 4. Blue Ball
    ball = balls[3]
    # View might be blocked by wardrobe. Best viewpoint is far left.
    view_pos = Point(joe_view.x_min, joe_view.y, joe_view.z_eye)
    obstruction_wardrobe_body = Point(9.5, 8, -1) # x=9.5, y=8
    
    print(f"4. Checking {ball.color} Ball at ({ball.pos.x}, {ball.pos.y}, {ball.pos.z})")
    print(f"   To see to the NE corner, the view might be blocked by the wardrobe body at x={obstruction_wardrobe_body.x}, y=[4, 8].")
    print(f"   Best viewpoint is far left: ({view_pos.x}, {view_pos.y}, {view_pos.z}).")

    vx, vy = view_pos.x, view_pos.y
    bx, by = ball.pos.x, ball.pos.y
    ox, oy = obstruction_wardrobe_body.x, obstruction_wardrobe_body.y

    # Line of sight equation in x-y plane
    sight_y_at_obstacle = vy + (ox - vx) * (by - vy) / (bx - vx)
    
    print(f"   Equation: y = {vy} + ({ox} - {vx}) * ({by} - {vy}) / ({bx} - {vx})")
    print(f"   Calculated y-position of sightline at wardrobe front (x={ox}): y = {sight_y_at_obstacle:.2f} ft.")
    print(f"   The sightline y-position ({sight_y_at_obstacle:.2f} ft) is greater than the wardrobe's top edge (y={oy} ft).")
    print(f"   Conclusion: {ball.color} ball is VISIBLE.\n")
    visible_balls.append(ball.color)

    # 5. Purple Ball
    ball = balls[4]
    # Ball is inside wardrobe. View blocked by wardrobe frame? Check from best viewpoint.
    view_pos = Point(joe_view.x_max, joe_view.y, joe_view.z_eye)
    obstruction_wardrobe_opening = Point(9.5, 4, -1) # x=9.5, y_min=4
    
    print(f"5. Checking {ball.color} Ball at ({ball.pos.x}, {ball.pos.y}, {ball.pos.z})")
    print(f"   Ball is inside the wardrobe. View is through the opening at x={obstruction_wardrobe_opening.x}, from y=4 to y=8.")
    print(f"   Using most optimistic viewpoint: ({view_pos.x}, {view_pos.y}, {view_pos.z}).")

    vx, vy = view_pos.x, view_pos.y
    px, py = ball.pos.x, ball.pos.y
    ox, oy = obstruction_wardrobe_opening.x, obstruction_wardrobe_opening.y

    # Line of sight equation in x-y plane
    sight_y_at_obstacle = vy + (ox - vx) * (py - vy) / (px - vx)
    
    print(f"   Equation: y = {vy} + ({ox} - {vx}) * ({py} - {vy}) / ({px} - {vx})")
    print(f"   Calculated y-position of sightline at wardrobe opening (x={ox}): y = {sight_y_at_obstacle:.2f} ft.")
    print(f"   The sightline y-position ({sight_y_at_obstacle:.2f} ft) is lower than the opening's bottom edge (y={oy} ft).")
    print(f"   Conclusion: {ball.color} ball is NOT VISIBLE (blocked by wardrobe frame).\n")
    
    print("--- Summary ---")
    print("The balls Joe can see are:")
    for b in visible_balls:
        print(f"- {b}")

analyze_visibility()
<<<Red, Yellow, Blue>>>