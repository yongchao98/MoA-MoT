import math

class Vec3:
    """A simple class for a 3D vector or point to make geometric operations cleaner."""
    def __init__(self, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __sub__(self, other):
        """Vector subtraction."""
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __repr__(self):
        """String representation for debugging."""
        return f"Vec3({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"

class AABB:
    """Represents an Axis-Aligned Bounding Box for use as an obstruction."""
    def __init__(self, name, min_corner, max_corner):
        self.name = name
        self.min = min_corner
        self.max = max_corner

    def intersects(self, ray_origin, ray_dir):
        """
        Checks for an intersection between a line segment and the box using the slab test.
        The line segment is defined by `ray_origin` and `ray_origin + ray_dir`.
        Returns True if the segment intersects the box, False otherwise.
        """
        t_min = 0.0
        t_max = 1.0  # We only care about intersections along the segment

        for axis in ['x', 'y', 'z']:
            origin_i = getattr(ray_origin, axis)
            dir_i = getattr(ray_dir, axis)
            min_i = getattr(self.min, axis)
            max_i = getattr(self.max, axis)

            if abs(dir_i) < 1e-6:  # Ray is parallel to the slab for this axis
                if origin_i < min_i or origin_i > max_i:
                    return False  # Parallel and outside the slab, so no intersection
            else:
                t1 = (min_i - origin_i) / dir_i
                t2 = (max_i - origin_i) / dir_i

                if t1 > t2:
                    t1, t2 = t2, t1  # Ensure t1 is the 'near' intersection

                t_min = max(t_min, t1)
                t_max = min(t_max, t2)

                if t_min >= t_max:
                    return False # Box is missed

        return True # Box is intersected by the line segment

def solve_visibility_problem():
    """
    This function models the room and objects in 3D to determine which balls Joe can see.
    It defines the coordinates of all items and performs line-of-sight checks against
    all potential obstructions for each ball from Joe's possible viewpoints.
    """
    # --- 1. Define Scene Geometry (all units in feet) ---

    # Joe's properties
    JOE_HEIGHT = 5.0
    JOE_LEAN_DEPTH = 3.0 / 12.0  # 3 inches converted to feet

    # Joe's possible viewpoints (extremes of the area he can occupy in the doorway)
    # Doorway is 3ft wide, centered on the 12ft south wall. Center is x=6.
    # Doorway runs from x = 6 - 1.5 = 4.5 to x = 6 + 1.5 = 7.5.
    joe_viewpoints = [
        Vec3(4.5, JOE_LEAN_DEPTH, JOE_HEIGHT),  # Westernmost viewpoint
        Vec3(7.5, JOE_LEAN_DEPTH, JOE_HEIGHT),  # Easternmost viewpoint
        Vec3(6.0, JOE_LEAN_DEPTH, JOE_HEIGHT),  # Center viewpoint
    ]

    # Ball properties
    BALL_RADIUS = (6.0 / 12.0) / 2.0  # 6-inch diameter -> 0.25 ft radius

    # Ball locations (Name -> Center Position)
    balls = {
        "Red": Vec3(12.0 - BALL_RADIUS, 0.0 + BALL_RADIUS, 0.0 + BALL_RADIUS),
        "Blue": Vec3(12.0 - BALL_RADIUS, 12.0 - BALL_RADIUS, 2.5 + BALL_RADIUS),
        "Yellow": Vec3(0.0 + BALL_RADIUS, 12.0 - BALL_RADIUS, 0.0 + BALL_RADIUS),
        "Green": Vec3(0.0 + BALL_RADIUS, 0.0 + BALL_RADIUS, 7.0 + BALL_RADIUS),
        "Purple": Vec3(12.0 - BALL_RADIUS, 4.0 + BALL_RADIUS, 0.0 + BALL_RADIUS),
    }

    # Obstruction properties, modeled as Axis-Aligned Bounding Boxes (AABBs)
    DOOR_HEIGHT = 80.0 / 12.0  # ~6.67 ft
    WARDROBE_HEIGHT = 6.0
    # Give open doors a tiny thickness to work with the AABB intersection model
    DOOR_THICKNESS = 0.01

    obstructions = [
        # Room door, hinged at x=4.5, open 90 degrees into the room, 3ft wide
        AABB("Room Door", Vec3(4.5, 0, 0), Vec3(4.5 + DOOR_THICKNESS, 3.0, DOOR_HEIGHT)),
        # Bookshelf on west wall (1ft deep, 4ft wide/long, 7ft tall)
        AABB("Bookshelf", Vec3(0, 0, 0), Vec3(1.0, 4.0, 7.0)),
        # Wardrobe main body on east wall (2.5ft deep, 4ft wide, centered)
        AABB("Wardrobe Body", Vec3(12.0 - 2.5, 4.0, 0), Vec3(12.0, 8.0, WARDROBE_HEIGHT)),
        # Wardrobe's open south door (2ft wide, hinged at (9.5, 4))
        AABB("Wardrobe South Door", Vec3(9.5, 2.0, 0), Vec3(9.5 + DOOR_THICKNESS, 4.0, WARDROBE_HEIGHT)),
        # Wardrobe's open north door (2ft wide, hinged at (9.5, 8))
        AABB("Wardrobe North Door", Vec3(9.5, 8.0, 0), Vec3(9.5 + DOOR_THICKNESS, 10.0, WARDROBE_HEIGHT)),
    ]

    # --- 2. Perform Visibility Check ---
    visible_balls = []
    print("Checking visibility for each ball...\n")
    for ball_name, ball_pos in balls.items():
        is_visible = False
        # Check for at least one clear line of sight from Joe's possible positions
        for joe_pos in joe_viewpoints:
            line_of_sight_vector = ball_pos - joe_pos
            
            is_blocked = False
            for obs in obstructions:
                if obs.intersects(joe_pos, line_of_sight_vector):
                    is_blocked = True
                    break  # This viewpoint is blocked, try the next viewpoint
            
            if not is_blocked:
                is_visible = True
                break  # This ball is visible, no need to check other viewpoints
        
        if is_visible:
            print(f"Result: The {ball_name} ball is VISIBLE.")
            visible_balls.append(ball_name)
        else:
            print(f"Result: The {ball_name} ball is NOT VISIBLE.")

    # --- 3. Output the Final Result ---
    # The problem asks for an "equation", but this is a logic problem.
    # The "equation" is the logical deduction based on the geometry.
    # Final equation: Visible_Balls = All_Balls - Blocked_Balls
    print("\n--- FINAL CALCULATION ---")
    print("This is a geometry problem solved by checking lines of sight against obstacles.")
    print("The final answer is derived from the following logical steps for each ball:")
    print("Is Visible(ball) = EXISTS(joe_pos) SUCH THAT NOT EXISTS(obstacle) WHERE Line(joe_pos, ball) intersects Obstacle")
    print("\nBased on this, the visible balls are determined.")
    
    if not visible_balls:
        final_answer_str = "Joe cannot see any of the balls."
    elif len(visible_balls) == 1:
        final_answer_str = f"Joe can see the {visible_balls[0]} ball."
    else:
        # Format output as "Ball1 and Ball2"
        visible_balls.sort()
        final_answer_str = " and ".join(visible_balls)
        final_answer_str = f"Joe can see the {final_answer_str} balls."

    print(f"\n{final_answer_str}")

    # Final answer in the required format
    final_answer_tag = " and ".join(sorted(visible_balls))
    print(f"\n<<<{final_answer_tag}>>>")

solve_visibility_problem()