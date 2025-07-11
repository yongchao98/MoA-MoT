import math

class Point:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z
    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)
    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y, self.z + other.z)
    def __mul__(self, scalar):
        return Point(self.x * scalar, self.y * scalar, self.z * scalar)
    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"

# Represents an Axis-Aligned Bounding Box
class AABB:
    def __init__(self, p_min, p_max):
        self.min = p_min
        self.max = p_max
    def __repr__(self):
        return f"AABB(min={self.min}, max={self.max})"

def check_line_aabb_intersection(p1, p2, box):
    """
    Checks if the line segment from p1 to p2 intersects with the AABB 'box'.
    Using the slab method.
    """
    direction = p2 - p1
    dir_inv = Point(
        1.0 / direction.x if direction.x != 0 else float('inf'),
        1.0 / direction.y if direction.y != 0 else float('inf'),
        1.0 / direction.z if direction.z != 0 else float('inf')
    )

    t_min = (box.min.x - p1.x) * dir_inv.x
    t_max = (box.max.x - p1.x) * dir_inv.x
    if t_min > t_max: t_min, t_max = t_max, t_min

    ty_min = (box.min.y - p1.y) * dir_inv.y
    ty_max = (box.max.y - p1.y) * dir_inv.y
    if ty_min > ty_max: ty_min, ty_max = ty_max, ty_min

    if (t_min > ty_max) or (ty_min > t_max):
        return False, "No intersection on XY plane", None
    
    t_min = max(t_min, ty_min)
    t_max = min(t_max, ty_max)

    tz_min = (box.min.z - p1.z) * dir_inv.z
    tz_max = (box.max.z - p1.z) * dir_inv.z
    if tz_min > tz_max: tz_min, tz_max = tz_max, tz_min

    if (t_min > tz_max) or (tz_min > t_max):
        return False, "No intersection on Z axis", None
        
    t_min = max(t_min, tz_min)
    t_max = min(t_max, tz_max)

    # Intersection point is p1 + direction * t_min
    intersection_point = p1 + direction * t_min
    
    # Check if intersection is within the segment [p1, p2] (i.e., t is between 0 and 1)
    if 0 <= t_min <= 1:
        return True, "Intersection found", intersection_point
        
    return False, "Intersection is outside the line segment", None


def analyze_visibility():
    """
    Main function to analyze the visibility of each ball.
    """
    # --- Room and Joe's Definition ---
    # Joe's eye position: he can move in the doorway [4.5, 7.5] and lean in 3 inches (0.25 ft)
    # Eye height is 5' - 3" = 4.75'
    joe_eye_height = 4.75
    # We check from the most advantageous points
    # V_left is better for seeing things on the west side, V_right for the east side
    vantage_points = {
        'left': Point(4.5, 0.25, joe_eye_height),
        'right': Point(7.5, 0.25, joe_eye_height),
    }

    # --- Obstacle Definitions (as AABBs) ---
    wardrobe_height = 8 # Assumed height
    door_height = 7 # Assumed height
    obstacles = {
        'bookshelf': AABB(Point(0, 0, 0), Point(1, 4, 7)),
        'wardrobe_main': AABB(Point(9.5, 4, 0), Point(12, 8, wardrobe_height)),
        'wardrobe_door_S': AABB(Point(7.5, 4, 0), Point(9.5, 4.01, wardrobe_height)), # South door (thin)
        # Hinge on right side (x=7.5) to see Yellow ball. Door swings to x=7.5, 0<=y<=3
        'door': AABB(Point(7.49, 0, 0), Point(7.5, 3, door_height)),
    }

    # --- Ball Definitions (Targets) ---
    ball_radius = 0.25 # 6-inch diameter / 2
    balls = {
        'Red': Point(12 - ball_radius, 0 + ball_radius, ball_radius),
        'Blue': Point(12 - ball_radius, 12 - ball_radius, 2.5 + ball_radius), # on 2.5ft table
        'Yellow': Point(0 + ball_radius, 12 - ball_radius, ball_radius),
        'Green': Point(0 + ball_radius, 0 + ball_radius, 7 + ball_radius), # on 7ft shelf
        'Purple': Point(12 - ball_radius, 4 + ball_radius, 0.5 + ball_radius), # inside wardrobe
    }
    
    visible_balls = []

    print("--- Starting Visibility Analysis ---")
    
    for ball_name, ball_pos in balls.items():
        print(f"\n--- Checking {ball_name} Ball at {ball_pos} ---")
        is_visible = False
        
        # Determine best vantage point to check from
        vantage_key = 'left' if ball_pos.x < 6 else 'right'
        eye_pos = vantage_points[vantage_key]
        
        print(f"Checking from vantage point {vantage_key}: {eye_pos}")

        is_blocked = False
        for obs_name, obstacle in obstacles.items():
            # A ball can't be blocked by an object it's on top of if we're looking from above.
            # However, for the green ball, the eye is below the ball, so the shelf itself can block.
            
            # For Purple, it's inside the wardrobe, so let's only check against the main body
            if ball_name == 'Purple' and obs_name != 'wardrobe_main':
                continue
            
            intersect, reason, point = check_line_aabb_intersection(eye_pos, ball_pos, obstacle)
            if intersect:
                is_blocked = True
                print(f"  [BLOCKED] by '{obs_name}'.")
                print(f"  Equation: P = Eye + t * (Ball - Eye)")
                print(f"  Eye = {eye_pos}, Ball = {ball_pos}")
                direction = ball_pos - eye_pos
                print(f"  P = {eye_pos} + t * {direction}")
                print(f"  Line of sight intersects obstacle at point {point} which is on the segment.")
                break
        
        if not is_blocked:
            is_visible = True
            
        if is_visible:
            print(f"  [VISIBLE] No direct obstructions found from this vantage point.")
            visible_balls.append(ball_name)
        else:
            print(f"  [NOT VISIBLE] View is obstructed.")
    
    print("\n--- Final Result ---")
    print(f"The balls Joe can see are: {', '.join(visible_balls) if visible_balls else 'None'}.")
    
    return visible_balls

visible_balls = analyze_visibility()
# The format requested by the prompt for the final answer.
# Printing the result in this format.
final_answer = ", ".join(sorted(visible_balls))
print(f"\n<<<{final_answer}>>>")
