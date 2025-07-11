import math

def solve_visibility():
    """
    Calculates which balls Joe can see based on the room's layout and obstructions.
    The room is a 12x12 space with the SW corner at (0,0).
    """
    
    # --- Define room elements and their coordinates ---
    
    # Joe's possible eye positions (x_range, y_range, z_height)
    joe_pos = {'x_min': 4.5, 'x_max': 7.5, 'y_max': 0.25, 'z': 4.75}
    
    # Obstructions
    # Bookshelf: x=[0,1], y=[0,4], z=7
    bookshelf = {'x_start': 0, 'x_end': 1, 'y_start': 0, 'y_end': 4, 'z': 7}
    # Wardrobe: x=[9.5,12], y=[4,8]
    wardrobe = {'x_start': 9.5, 'x_end': 12, 'y_start': 4, 'y_end': 8, 'z': 6} # Assume 6ft tall
    # Open Door: 3ft wide, so 3ft high when perpendicular
    door_height = 3
    
    # Ball positions (x, y, z)
    balls = {
        'red':    {'pos': (11.75, 0.25, 0.25), 'visible': False, 'reason': ''},
        'blue':   {'pos': (11.75, 11.75, 3.25), 'visible': False, 'reason': ''},
        'yellow': {'pos': (0.25, 11.75, 0.25), 'visible': False, 'reason': ''},
        'green':  {'pos': (0.25, 0.25, 7.25), 'visible': False, 'reason': ''},
        'purple': {'pos': (11.75, 4.25, 0.25), 'visible': False, 'reason': ''}
    }
    
    visible_balls = []

    print("Analyzing visibility for each ball:\n")

    # --- Analysis for each ball ---

    # 1. YELLOW BALL (Northwest Corner)
    # To see around the bookshelf (at x=1), Joe stands as far right as possible.
    # We assume the door is hinged at x=7.5, so it doesn't block the view to the left.
    joe_viewpoint = (joe_pos['x_max'], joe_pos['y_max'])
    ball_pos = balls['yellow']['pos']
    bookshelf_edge_x = bookshelf['x_end']
    slope = (ball_pos[1] - joe_viewpoint[1]) / (ball_pos[0] - joe_viewpoint[0])
    y_intersect = joe_viewpoint[1] + slope * (bookshelf_edge_x - joe_viewpoint[0])
    
    balls['yellow']['reason'] = (
        f"To see the Yellow ball at {ball_pos[:2]}, Joe stands at the right of the doorway ({joe_viewpoint}).\n"
        f"The line of sight must clear the bookshelf, which ends at y = {bookshelf['y_end']}.\n"
        f"The equation for the line of sight is y - {joe_viewpoint[1]} = (({ball_pos[1]} - {joe_viewpoint[1]}) / ({ball_pos[0]} - {joe_viewpoint[0]})) * (x - {joe_viewpoint[0]}).\n"
        f"At the bookshelf's edge (x = {bookshelf_edge_x}), the line's height is y = {y_intersect:.2f}.\n"
        f"Since {y_intersect:.2f} > {bookshelf['y_end']}, the line of sight passes the bookshelf."
    )
    if y_intersect > bookshelf['y_end']:
        balls['yellow']['visible'] = True
        
    # 2. GREEN BALL (On Bookshelf)
    # Check if the line of sight in the Z-dimension clears the bookshelf's top edge.
    joe_viewpoint_xz = (joe_pos['x_max'], joe_pos['z'])
    ball_pos_xz = (balls['green']['pos'][0], balls['green']['pos'][2])
    bookshelf_edge_xz = (bookshelf['x_end'], bookshelf['z'])
    slope_xz = (ball_pos_xz[1] - joe_viewpoint_xz[1]) / (ball_pos_xz[0] - joe_viewpoint_xz[0])
    z_intersect = joe_viewpoint_xz[1] + slope_xz * (bookshelf_edge_xz[0] - joe_viewpoint_xz[0])

    balls['green']['reason'] = (
        f"The Green ball is on the 7-foot-tall bookshelf at z = {balls['green']['pos'][2]}.\n"
        f"Joe's eye height is z = {joe_pos['z']}. He views from x = {joe_pos['x_max']}.\n"
        f"The line of sight must pass over the front edge of the bookshelf at (x={bookshelf_edge_xz[0]}, z={bookshelf_edge_xz[1]}).\n"
        f"The equation for the line of sight in the X-Z plane is z - {joe_viewpoint_xz[1]} = (({ball_pos_xz[1]} - {joe_viewpoint_xz[1]}) / ({ball_pos_xz[0]} - {joe_viewpoint_xz[0]})) * (x - {joe_viewpoint_xz[0]}).\n"
        f"At the bookshelf's edge (x = {bookshelf_edge_xz[0]}), the line's height is z = {z_intersect:.2f}.\n"
        f"Since {z_intersect:.2f} < {bookshelf_edge_xz[1]}, the line of sight is blocked by the top of the bookshelf."
    )
    if z_intersect > bookshelf_edge_xz[1]:
        balls['green']['visible'] = True

    # 3. RED BALL (Southeast Corner)
    # To see the Red ball, we assume the door is hinged at x=4.5, clearing the view to the right.
    # The line of sight is low, below the wardrobe.
    ball_y = balls['red']['pos'][1]
    wardrobe_y = wardrobe['y_start']
    balls['red']['reason'] = (
        f"Assuming the door hinge allows a clear view to the right, Joe can see the Red ball at y = {ball_y}.\n"
        f"The wardrobe, the only potential obstruction, starts at y = {wardrobe_y}.\n"
        f"The line of sight y = {ball_y} is below the wardrobe (y = {wardrobe_y})."
    )
    if ball_y < wardrobe_y:
        balls['red']['visible'] = True
        
    # 4. BLUE BALL (Northeast Corner)
    # To see around the wardrobe (at x=9.5), Joe stands as far left as possible.
    # We assume the door is hinged at x=7.5 to see to the right, but the view must also pass over the door.
    # First, check if view clears the wardrobe.
    joe_viewpoint = (joe_pos['x_min'], joe_pos['y_max'])
    ball_pos = balls['blue']['pos']
    wardrobe_edge_x = wardrobe['x_start']
    slope = (ball_pos[1] - joe_viewpoint[1]) / (ball_pos[0] - joe_viewpoint[0])
    y_intersect = joe_viewpoint[1] + slope * (wardrobe_edge_x - joe_viewpoint[0])
    
    balls['blue']['reason'] = (
        f"To see the Blue ball at {ball_pos[:2]}, Joe stands at the left of the doorway ({joe_viewpoint}).\n"
        f"The line of sight must clear the wardrobe, which is at y <= {wardrobe['y_end']} for x > {wardrobe_edge_x}.\n"
        f"The equation for the line of sight is y - {joe_viewpoint[1]} = (({ball_pos[1]} - {joe_viewpoint[1]}) / ({ball_pos[0]} - {joe_viewpoint[0]})) * (x - {joe_viewpoint[0]}).\n"
        f"At the wardrobe's edge (x = {wardrobe_edge_x}), the line's y-coordinate is {y_intersect:.2f}.\n"
        f"Since {y_intersect:.2f} > {wardrobe['y_end']}, the line of sight passes the wardrobe body.\n"
    )
    # Second, check if the view clears the open door (hinged at x=7.5).
    door_plane_x = joe_pos['x_max']
    y_intersect_door = joe_viewpoint[1] + slope * (door_plane_x - joe_viewpoint[0])
    balls['blue']['reason'] += (
        f"The open door (hinged at {door_plane_x}) is {door_height} feet high. The line of sight crosses the door's plane at y={y_intersect_door:.2f}, which is above {door_height} feet."
    )
    if y_intersect > wardrobe['y_end'] and y_intersect_door > door_height:
        balls['blue']['visible'] = True

    # 5. PURPLE BALL (In Wardrobe)
    # Check if the line of sight can enter the wardrobe opening (y > 4 at x=9.5).
    joe_viewpoint = (joe_pos['x_min'], joe_pos['y_max'])
    ball_pos = balls['purple']['pos']
    wardrobe_edge_x = wardrobe['x_start']
    wardrobe_opening_y = wardrobe['y_start']
    slope = (ball_pos[1] - joe_viewpoint[1]) / (ball_pos[0] - joe_viewpoint[0])
    y_intersect = joe_viewpoint[1] + slope * (wardrobe_edge_x - joe_viewpoint[0])

    balls['purple']['reason'] = (
        f"The Purple ball is inside the wardrobe at {ball_pos[:2]}.\n"
        f"The wardrobe opening starts at y = {wardrobe_opening_y} at x = {wardrobe_edge_x}.\n"
        f"From the best vantage point ({joe_viewpoint}), the equation for the line of sight is y - {joe_viewpoint[1]} = (({ball_pos[1]} - {joe_viewpoint[1]}) / ({ball_pos[0]} - {joe_viewpoint[0]})) * (x - {joe_viewpoint[0]}).\n"
        f"This line intersects the wardrobe's front (x = {wardrobe_edge_x}) at y = {y_intersect:.2f}.\n"
        f"Since {y_intersect:.2f} < {wardrobe_opening_y}, the line of sight hits the solid part of the wardrobe."
    )
    if y_intersect > wardrobe_opening_y:
        balls['purple']['visible'] = True

    # --- Print final results ---
    print("-" * 20)
    for ball_name, data in balls.items():
        visibility = "Visible" if data['visible'] else "Not Visible"
        print(f"Ball: {ball_name.capitalize()}")
        print(f"Result: {visibility}")
        print(f"Reasoning: {data['reason']}\n")
        if data['visible']:
            visible_balls.append(ball_name.capitalize())
    
    print("-" * 20)
    print("Conclusion:")
    print(f"Joe can see the following balls: {', '.join(visible_balls)}.")
    
solve_visibility()
<<<Red, Blue, Yellow>>>