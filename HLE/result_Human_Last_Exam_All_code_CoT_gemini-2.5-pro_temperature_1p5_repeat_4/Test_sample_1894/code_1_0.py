import math

class Point3D:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z
    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"

class Box3D:
    """An axis-aligned 3D box obstacle."""
    def __init__(self, name, x_range, y_range, z_range):
        self.name = name
        self.min_p = Point3D(x_range[0], y_range[0], z_range[0])
        self.max_p = Point3D(x_range[1], y_range[1], z_range[1])

    def intersects(self, p1, p2):
        """
        Checks if the line segment from p1 to p2 intersects the box.
        Uses the slab method for line-box intersection.
        """
        # Direction vector of the line
        dir_vec = Point3D(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
        
        # Handle division by zero for axis-aligned lines
        epsilon = 1e-6
        dir_vec.x = epsilon if abs(dir_vec.x) < epsilon else dir_vec.x
        dir_vec.y = epsilon if abs(dir_vec.y) < epsilon else dir_vec.y
        dir_vec.z = epsilon if abs(dir_vec.z) < epsilon else dir_vec.z
        
        # Calculate intersection t-values for each pair of slab planes
        t_near = -math.inf
        t_far = math.inf

        # X-slabs
        t1 = (self.min_p.x - p1.x) / dir_vec.x
        t2 = (self.max_p.x - p1.x) / dir_vec.x
        if t1 > t2: t1, t2 = t2, t1
        t_near = max(t_near, t1)
        t_far = min(t_far, t2)

        # Y-slabs
        t1 = (self.min_p.y - p1.y) / dir_vec.y
        t2 = (self.max_p.y - p1.y) / dir_vec.y
        if t1 > t2: t1, t2 = t2, t1
        t_near = max(t_near, t1)
        t_far = min(t_far, t2)

        # Z-slabs
        t1 = (self.min_p.z - p1.z) / dir_vec.z
        t2 = (self.max_p.z - p1.z) / dir_vec.z
        if t1 > t2: t1, t2 = t2, t1
        t_near = max(t_near, t1)
        t_far = min(t_far, t2)
        
        # A true intersection occurs if t_near < t_far
        # and the intersection interval [t_near, t_far] overlaps with the segment interval [0, 1]
        if t_near > t_far or t_far < 0 or t_near > 1:
            return None # No intersection
            
        intersection_t = t_near
        intersection_point = Point3D(
            p1.x + intersection_t * dir_vec.x,
            p1.y + intersection_t * dir_vec.y,
            p1.z + intersection_t * dir_vec.z
        )
        
        return intersection_point


def solve():
    """
    Solves the visibility problem for Joe in the room.
    """
    # Define obstacles in the room
    # Note: Door/shelf heights are assumed based on typical dimensions
    obstacles = [
        Box3D("Bookshelf", x_range=(0, 1), y_range=(0, 4), z_range=(0, 7)),
        Box3D("Wardrobe", x_range=(9.5, 12), y_range=(4, 8), z_range=(0, 6)),
        Box3D("Wardrobe Left Door", x_range=(7.5, 9.5), y_range=(4, 4.1), z_range=(0, 6)),
        Box3D("Wardrobe Right Door", x_range=(7.5, 9.5), y_range=(7.9, 8), z_range=(0, 6)),
        Box3D("Room Door", x_range=(4.5, 4.6), y_range=(0, 3), z_range=(0, 7))
    ]

    # Define ball locations
    balls = {
        "Red": Point3D(11.9, 0.1, 0.3),  # On floor in SE corner
        "Blue": Point3D(11.9, 11.9, 2.5), # On table in NE corner
        "Yellow": Point3D(0.1, 11.9, 0.3),# On floor in NW corner
        "Green": Point3D(1.0, 4.0, 7.3),  # On top of the bookshelf
        "Purple": Point3D(11.9, 4.2, 0.3) # Inside wardrobe
    }

    # Define Joe's possible eye positions (test key corners of viewing area)
    # y = -0.25 corresponds to leaning in 3 inches
    joe_view_points = [
        Point3D(4.5, -0.25, 4.75), # Farthest left and forward
        Point3D(7.5, -0.25, 4.75), # Farthest right and forward
    ]

    visible_balls = []

    print("Analyzing visibility for each ball:\n")

    for color, ball_pos in balls.items():
        is_visible = False
        print(f"--- Checking {color} Ball at {ball_pos} ---")
        
        for i, joe_pos in enumerate(joe_view_points):
            blocker = None
            for obs in obstacles:
                intersection_point = obs.intersects(joe_pos, ball_pos)
                if intersection_point:
                    # Final check: is the intersection point actually between Joe and the ball?
                    dist_joe_ball = math.dist((joe_pos.x, joe_pos.y, joe_pos.z), (ball_pos.x, ball_pos.y, ball_pos.z))
                    dist_joe_intersect = math.dist((joe_pos.x, joe_pos.y, joe_pos.z), (intersection_point.x, intersection_point.y, intersection_point.z))
                    if dist_joe_intersect < dist_joe_ball - 1e-6: # Add tolerance
                         blocker = (obs, intersection_point)
                         break # Blocked by this obstacle
            
            if not blocker:
                print(f"  Line of sight from Joe at {joe_pos} is CLEAR.")
                is_visible = True
                break # Found a clear line of sight, no need to check other view points
            else:
                 obs, p = blocker
                 print(f"  Line of sight from Joe at {joe_pos} is BLOCKED by {obs.name}.")
                 equation_msg = (
                    f"    Joe's Position (P1): ({joe_pos.x}, {joe_pos.y}, {joe_pos.z}), "
                    f"Ball's Position (P2): ({ball_pos.x}, {ball_pos.y}, {ball_pos.z})\n"
                    f"    The line between P1 and P2 intersects the obstacle at approx. ({p.x:.2f}, {p.y:.2f}, {p.z:.2f})."
                 )
                 print(equation_msg)

        if is_visible:
            print(f"Conclusion: The {color} ball is VISIBLE.\n")
            visible_balls.append(color)
        else:
            print(f"Conclusion: The {color} ball is NOT VISIBLE from any tested viewpoint.\n")

    print("--- Final Result ---")
    print(f"Joe can see the following balls: {', '.join(visible_balls)}.")


solve()
<<<Red, Yellow>>>