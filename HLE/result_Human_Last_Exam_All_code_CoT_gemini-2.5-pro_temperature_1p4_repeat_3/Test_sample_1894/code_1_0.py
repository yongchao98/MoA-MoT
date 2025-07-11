import numpy as np

def get_line_y_at_x(p1, p2, x):
    """Calculates the y-coordinate on the line defined by p1-p2 at a given x."""
    if p2[0] - p1[0] == 0: return p1[1] # Vertical line
    m = (p2[1] - p1[1]) / (p2[0] - p1[0])
    c = p1[1] - m * p1[0]
    return m * x + c

def get_line_x_at_y(p1, p2, y):
    """Calculates the x-coordinate on the line defined by p1-p2 at a given y."""
    if p2[1] - p1[1] == 0: return p1[0] # Horizontal line
    m_inv = (p2[0] - p1[0]) / (p2[1] - p1[1])
    c_inv = p1[0] - m_inv * p1[1]
    return m_inv * y + c_inv

# --- Define Points ---
# Joe's eye can be anywhere on the line segment p1-p2
joe_p1 = {'name': "left doorway edge", 'x': 4.5, 'y': 0.25, 'z': 5.0}
joe_p2 = {'name': "right doorway edge", 'x': 7.5, 'y': 0.25, 'z': 5.0}

balls = {
    'Red':    {'x': 11.75, 'y': 0.25,  'z': 0.25,  'visible': False, 'reason': ''},
    'Green':  {'x': 0.25,  'y': 0.25,  'z': 7.25,  'visible': False, 'reason': ''},
    'Yellow': {'x': 0.25,  'y': 11.75, 'z': 0.25,  'visible': True,  'reason': ''},
    'Blue':   {'x': 11.75, 'y': 11.75, 'z': 2.75,  'visible': False, 'reason': ''},
    'Purple': {'x': 11.75, 'y': 4.25,  'z': 0.25,  'visible': False, 'reason': ''}
}

print("Analyzing visibility for each ball...\n")

# --- Analysis ---

# 1. Green Ball
# Blocked by the bookshelf it is sitting on. Joe's eye (z=5) is below the shelf top (z=7).
# The line of sight in the XY plane also passes through the bookshelf's footprint.
y_at_shelf_front = get_line_y_at_x(joe_p1, balls['Green'], 1.0)
balls['Green']['reason'] = f"Blocked by the 7-foot-tall bookshelf it sits on. Joe's eye level is 5 feet. The line of sight from the {joe_p1['name']} ({joe_p1['x']}, {joe_p1['y']}) to the ball ({balls['Green']['x']}, {balls['Green']['y']}) has a y-coordinate of {y_at_shelf_front:.2f} at the bookshelf's front (x=1.0), which is within the shelf's footprint (y from 0 to 4)."
balls['Green']['visible'] = False

# 2. Red Ball
# The line of sight is at a constant y=0.25. This line must pass through the open door plane at x=7.5.
# The door occupies y from 0 to 3. Since y=0.25 is in this range, the door blocks the view.
balls['Red']['reason'] = f"Blocked by the open door. The line of sight from any point in the doorway to the ball is at y=0.25. This line must pass through the plane of the door at x=7.5. Since the y-coordinate 0.25 is within the door's range (y from 0 to 3), the door obstructs the view."
balls['Red']['visible'] = False

# 3. Yellow Ball
# Must check for obstruction by the bookshelf corner at (1, 4).
# Check from Joe's left-most view
y_at_shelf_corner = get_line_y_at_x(joe_p1, balls['Yellow'], 1.0)
if y_at_shelf_corner <= 4:
    balls['Yellow']['visible'] = False
    balls['Yellow']['reason'] = f"The view from the {joe_p1['name']} is blocked by the bookshelf, as the line of sight passes too close to its corner."
# Check from Joe's right-most view
y_at_shelf_corner_2 = get_line_y_at_x(joe_p2, balls['Yellow'], 1.0)
if y_at_shelf_corner_2 <= 4:
    balls['Yellow']['visible'] = False
    balls['Yellow']['reason'] = f"The view from the {joe_p2['name']} is blocked by the bookshelf, as the line of sight passes too close to its corner."
else:
    balls['Yellow']['reason'] = f"Clear line of sight. From the {joe_p1['name']}, the sightline passes the bookshelf's front (x=1.0) at y={y_at_shelf_corner:.2f}, well clear of its corner (y=4). From the {joe_p2['name']}, it passes at y={y_at_shelf_corner_2:.2f}, also clear."

# 4. Blue Ball
# Check for obstruction by the two open wardrobe doors.
# View from Joe's left side to the Blue ball.
x_at_wardrobe_door2 = get_line_x_at_y(joe_p1, balls['Blue'], 8.0)
if 7.5 <= x_at_wardrobe_door2 <= 9.5:
    balls['Blue']['visible'] = False
    balls['Blue']['reason'] = f"Blocked by the upper wardrobe door (at y=8.0). The line of sight from the {joe_p1['name']} intersects this door's plane at x={x_at_wardrobe_door2:.2f}, which is within the door's width (x from 7.5 to 9.5)."
else: # If not blocked from the left, check the right side
    x_at_wardrobe_door1 = get_line_x_at_y(joe_p2, balls['Blue'], 4.0)
    if 7.5 <= x_at_wardrobe_door1 <= 9.5:
        balls['Blue']['visible'] = False
        balls['Blue']['reason'] = f"Blocked by the lower wardrobe door (at y=4.0). The line of sight from the {joe_p2['name']} intersects this door's plane at x={x_at_wardrobe_door1:.2f}, which is within the door's width (x from 7.5 to 9.5)."
    else: # If still not blocked (which is unlikely), it might be visible
        balls['Blue']['visible'] = True
        balls['Blue']['reason'] = "Clear line of sight."


# 5. Purple Ball
# The ball is inside the wardrobe. The view is only possible through the opening at x=9.5, between y=4 and y=8.
# Find where the line of sight hits the wardrobe's front plane (x=9.5)
y_at_wardrobe_front_1 = get_line_y_at_x(joe_p1, balls['Purple'], 9.5)
y_at_wardrobe_front_2 = get_line_y_at_x(joe_p2, balls['Purple'], 9.5)
if not (4 <= y_at_wardrobe_front_1 <= 8) and not (4 <= y_at_wardrobe_front_2 <= 8):
     balls['Purple']['visible'] = False
     balls['Purple']['reason'] = f"Blocked by the front of the wardrobe. The line of sight from the {joe_p1['name']} reaches the wardrobe's front plane (x=9.5) at y={y_at_wardrobe_front_1:.2f}, which is below the opening (y from 4 to 8)."
else:
     balls['Purple']['visible'] = True # This case won't be reached
     balls['Purple']['reason'] = "Clear line of sight through the wardrobe opening."


# --- Final Output ---
visible_balls = []
for name, data in balls.items():
    status = "VISIBLE" if data['visible'] else "NOT VISIBLE"
    print(f"- {name} Ball: {status}")
    print(f"  Reason: {data['reason']}\n")
    if data['visible']:
        visible_balls.append(name)

print("Final Answer: The only ball Joe can see is:")
for ball in visible_balls:
    print(f">>>{ball}<<<")