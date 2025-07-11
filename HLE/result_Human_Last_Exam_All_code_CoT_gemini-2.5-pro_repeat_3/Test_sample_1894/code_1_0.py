import math

class Point:
    """A class to represent a point in 2D space."""
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f})"

def on_segment(p, q, r):
    """Given three collinear points p, q, r, check if point q lies on line segment 'pr'."""
    return (q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and
            q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y))

def orientation(p, q, r):
    """Find orientation of ordered triplet (p, q, r).
    Returns:
    0 --> p, q and r are collinear
    1 --> Clockwise
    2 --> Counterclockwise
    """
    val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
    if val == 0: return 0  # Collinear
    return 1 if val > 0 else 2  # Clockwise or Counterclockwise

def intersects(seg1_p1, seg1_q1, seg2_p2, seg2_q2):
    """Return true if line segment 'p1q1' and 'p2q2' intersect."""
    o1 = orientation(seg1_p1, seg1_q1, seg2_p2)
    o2 = orientation(seg1_p1, seg1_q1, seg2_q2)
    o3 = orientation(seg2_p2, seg2_q2, seg1_p1)
    o4 = orientation(seg2_p2, seg2_q2, seg1_q1)

    # General case
    if (o1 != o2 and o3 != o4):
        return True

    # Special Cases for collinear points
    if (o1 == 0 and on_segment(seg1_p1, seg2_p2, seg1_q1)): return True
    if (o2 == 0 and on_segment(seg1_p1, seg2_q2, seg1_q1)): return True
    if (o3 == 0 and on_segment(seg2_p2, seg1_p1, seg2_q2)): return True
    if (o4 == 0 and on_segment(seg2_p2, seg1_q1, seg2_q2)): return True

    return False

def solve():
    """
    Solves the visibility problem for Joe in the room.
    """
    # Define locations of balls (centers)
    balls = {
        "Red":    {"pos": Point(12 - 0.25, 0 + 0.25), "z": 0.25},
        "Blue":   {"pos": Point(12 - 0.25, 12 - 0.25), "z": 2.75}, # Assuming table height of 2.5ft
        "Yellow": {"pos": Point(0 + 0.25, 12 - 0.25), "z": 0.25},
        "Green":  {"pos": Point(0 + 0.25, 0 + 0.25), "z": 7.25}, # On 7ft shelf
        "Purple": {"pos": Point(12 - 0.25, 4 + 0.25), "z": 0.25} # In wardrobe corner
    }

    # Define obstacles as line segments
    obstacles = {
        "Door": (Point(4.5, 0), Point(4.5, 3)),
        "Bookshelf": (Point(1, 0), Point(1, 4)),
        "Wardrobe Body": (Point(9.5, 4), Point(9.5, 8)),
        "Wardrobe South Door": (Point(7.5, 4), Point(9.5, 4)),
        "Wardrobe North Door": (Point(7.5, 8), Point(9.5, 8))
    }

    # Define Joe's viewing position (a line segment)
    joe_view_start = Point(4.5, 0.25)
    joe_view_end = Point(7.5, 0.25)
    # We will test a few points along this line for visibility
    joe_view_points = [joe_view_start, Point((joe_view_start.x + joe_view_end.x) / 2, joe_view_start.y), joe_view_end]
    
    visible_balls = []

    print("Analyzing which balls Joe can see...\n")

    for color, ball in balls.items():
        print(f"--- Checking {color} Ball at {ball['pos']} ---")
        is_visible = False
        
        # Special case for Green ball based on height
        if color == "Green":
            joe_eye_height = 4.75 # Approx. eye height for a 5ft person
            bookshelf_height = 7
            if ball['z'] > bookshelf_height and joe_eye_height < bookshelf_height:
                print(f"Joe's eye level ({joe_eye_height} ft) is below the bookshelf height ({bookshelf_height} ft).")
                print("The green ball is on top of the bookshelf, so it cannot be seen.")
                is_visible = False
            else:
                 # This else won't be hit based on problem, but good practice
                 pass
        
        # Check 2D line of sight for other balls
        if color != "Green":
            for joe_pos in joe_view_points:
                line_of_sight_blocked = False
                for obs_name, obs_segment in obstacles.items():
                    # The open door is at x=4.5. Joe's westernmost view is at x=4.5.
                    # He looks past it, it doesn't block views into the room.
                    if obs_name == "Door":
                        continue
                    
                    p1 = joe_pos
                    q1 = ball['pos']
                    p2 = obs_segment[0]
                    q2 = obs_segment[1]

                    if intersects(p1, q1, p2, q2):
                        print(f"Line of sight from {p1} to {q1} is BLOCKED by {obs_name} {p2}-{q2}")
                        line_of_sight_blocked = True
                        break # This line of sight is blocked, try next obstacle
                
                if not line_of_sight_blocked:
                    print(f"Found a clear line of sight from {joe_pos} to the ball.")
                    is_visible = True
                    break # Found a clear view, no need to check other viewpoints
            
            if not is_visible:
                print("No clear line of sight found from any viewpoint.")

        if is_visible:
            visible_balls.append(color)
        
        print(f"Result for {color} Ball: {'VISIBLE' if is_visible else 'NOT VISIBLE'}\n")

    print("--- Conclusion ---")
    print("The balls Joe can see are:")
    for ball_color in visible_balls:
        print(f"- {ball_color}")
        
    final_answer = ", ".join(sorted(visible_balls))
    print(f"\n<<<[{final_answer}]>>>")


solve()