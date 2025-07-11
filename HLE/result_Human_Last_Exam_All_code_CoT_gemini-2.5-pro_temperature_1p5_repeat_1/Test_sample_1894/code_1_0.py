import math

# A function to check the orientation of an ordered triplet (p, q, r)
# 0 --> p, q and r are collinear
# 1 --> Clockwise
# 2 --> Counterclockwise
def orientation(p, q, r):
    val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
    if val == 0: return 0  # Collinear
    return 1 if val > 0 else 2  # Clockwise or Counterclockwise

# Given three collinear points p, q, r, the function checks if point q lies on line segment 'pr'
def on_segment(p, q, r):
    return (q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and
            q[1] <= max(p[1], r[1]) and q[1] >= min(p[1], r[1]))

# The main function to find if line segment 'p1q1' and 'p2q2' intersect.
def intersects(seg1, seg2):
    p1, q1 = seg1
    p2, q2 = seg2

    # Find the four orientations needed for general and special cases
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case
    if (o1 != o2 and o3 != o4):
        return True

    # Special Cases for collinear points
    if (o1 == 0 and on_segment(p1, p2, q1)): return True
    if (o2 == 0 and on_segment(p1, q2, q1)): return True
    if (o3 == 0 and on_segment(p2, p1, q2)): return True
    if (o4 == 0 and on_segment(p2, q1, q2)): return True

    return False

def solve():
    # --- Define coordinates ---
    # Joe's most advantageous viewpoint for a wide view
    joe_pos = (4.5, 0.25)

    # Ball positions (x, y)
    balls = {
        'red': {'pos': (11.75, 0.25), 'desc': 'in the southeast corner'},
        'blue': {'pos': (11.75, 11.75), 'desc': 'in the far northeast corner'},
        'yellow': {'pos': (0.25, 11.75), 'desc': 'in the northwest corner'},
        'green': {'pos': (0.25, 0.25), 'desc': 'in the southwest corner on top of the shelf'},
        'purple': {'pos': (11.75, 4.25), 'desc': 'in the southeast corner of the wardrobe'}
    }

    # Obstruction line segments
    obstructions = {
        "Bookshelf body (solid)": {'type': 'area', 'x': (0, 1), 'y': (0, 4)},
        "Wardrobe setup (solid)": {'type': 'area', 'x': (7.5, 12), 'y': (4, 8)},
        "Open Door": {'type': 'line', 'seg': ((7.5, 0), (7.5, 3))},
        "Wardrobe South Door": {'type': 'line', 'seg': ((7.5, 4), (9.5, 4))},
        "Wardrobe North Door": {'type': 'line', 'seg': ((7.5, 8), (9.5, 8))},
        "Wardrobe Back": {'type': 'line', 'seg': ((9.5, 4), (9.5, 8))}
    }

    visible_balls = []

    print(f"Analyzing visibility from Joe's position at {joe_pos}\n")

    for color, ball in balls.items():
        ball_pos = ball['pos']
        is_visible = True
        blocker_name = None
        blocker_details = ""

        print(f"--- Checking {color.capitalize()} ball {ball['desc']} at {ball_pos} ---")

        # Initial placement checks
        if obstructions["Bookshelf body (solid)"]['x'][0] <= ball_pos[0] <= obstructions["Bookshelf body (solid)"]['x'][1] and \
           obstructions["Bookshelf body (solid)"]['y'][0] <= ball_pos[1] <= obstructions["Bookshelf body (solid)"]['y'][1]:
            print(f"Result: The Green ball is on the bookshelf which is between Joe and the ball's location.")
            is_visible = False
            blocker_name = "Bookshelf body"

        if is_visible and ball_pos[0] > 9.5 and 4 <= ball_pos[1] <= 8:
             print(f"Result: The Purple ball is located inside the footprint of the wardrobe.")
             is_visible = False
             blocker_name = "Wardrobe Body"

        if is_visible:
            line_of_sight = (joe_pos, ball_pos)
            for name, obs in obstructions.items():
                if obs['type'] == 'line':
                    if intersects(line_of_sight, obs['seg']):
                        is_visible = False
                        blocker_name = name
                        p1, p2 = obs['seg']
                        
                        # Add equation details
                        jx, jy = joe_pos
                        bx, by = ball_pos
                        # Line equation: (y - jy) = m * (x - jx) -> m = (by - jy) / (bx - jx)
                        if (bx - jx) != 0:
                            m = (by - jy) / (bx - jx)
                            # Find intersection point
                            if p1[0] == p2[0]: # Vertical obstruction line
                                intersect_x = p1[0]
                                intersect_y = m * (intersect_x - jx) + jy
                                blocker_details = f"Line of sight from {joe_pos} to {ball_pos} intersects the vertical line x = {intersect_x:.2f} at y = {intersect_y:.2f}."
                            else: # Horizontal obstruction line
                                intersect_y = p1[1]
                                if m != 0:
                                    intersect_x = (intersect_y - jy) / m + jx
                                    blocker_details = f"Line of sight from {joe_pos} to {ball_pos} intersects the horizontal line y = {intersect_y:.2f} at x = {intersect_x:.2f}."
                                else: # Both lines horizontal
                                    blocker_details = f"Line of sight is horizontal (y={jy}) and obstructed by horizontal line (y={intersect_y})."
                        else: # vertical line of sight
                            blocker_details = f"Line of sight is vertical (x={jx})."

                        print(f"Result: View is blocked by the '{blocker_name}' from {p1} to {p2}.")
                        print(f"Calculation: {blocker_details}")
                        break
        
        if is_visible:
            visible_balls.append(color.capitalize())
            print("Result: No direct obstructions found.")

        print(f"Conclusion: {color.capitalize()} ball is {'VISIBLE' if is_visible else 'NOT VISIBLE'}.\n")
    
    print("--- FINAL SUMMARY ---")
    if not visible_balls:
        print("Joe cannot see any of the balls.")
    else:
        print(f"Joe can see the following ball(s): {', '.join(visible_balls)}")


solve()
<<<Yellow>>>