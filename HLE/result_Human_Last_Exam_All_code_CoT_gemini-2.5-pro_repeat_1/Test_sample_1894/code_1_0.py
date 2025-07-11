import math

# Define a simple Point class for clarity
class Point:
    """A helper class to represent a 3D point and perform vector operations."""
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"
    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)
    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y, self.z + other.z)
    def __mul__(self, scalar):
        return Point(self.x * scalar, self.y * scalar, self.z * scalar)

def solve():
    """
    Determines which balls Joe can see by checking lines of sight.
    """
    print("Analyzing visibility for Joe...\n")

    # --- Define Scene Geometry ---
    # Joe's eye level (5ft person) and position
    JOE_EYE_Z = 4.75
    V_LEFT = Point(4.5, -0.25, JOE_EYE_Z)
    V_RIGHT = Point(7.5, -0.25, JOE_EYE_Z)

    # Ball positions (center of a 0.5ft diameter ball, radius = 0.25 ft)
    BALLS = {
        "Red":    Point(12 - 0.25, 0 + 0.25, 0.25),            # SE corner, on floor
        "Blue":   Point(12 - 0.25, 12 - 0.25, 2.5 + 0.25),     # NE corner, on 2.5ft table
        "Yellow": Point(0 + 0.25, 12 - 0.25, 0.25),            # NW corner, on floor
        "Green":  Point(0 + 0.25, 0 + 0.25, 7 + 0.25),         # SW corner, on 7ft shelf
        "Purple": Point(12 - 0.25, 4 + 0.25, 0.25)             # In SE corner of wardrobe
    }
    
    visible_balls = []

    # --- Red Ball (SE corner) ---
    print("Red ball (SE corner):")
    print(f"  - The Red ball is at {BALLS['Red']}.")
    print("  - The line of sight is towards the southeast corner. The wardrobe, the main potential obstacle,")
    print("    is positioned against the east wall starting at y=4. The ball is at y=0.25.")
    print("  - The view is clearly unobstructed.")
    print("  - Visibility: Visible\n")
    visible_balls.append("Red")

    # --- Purple Ball (in wardrobe) ---
    print("Purple ball (in wardrobe):")
    target = BALLS['Purple']
    print(f"  - The Purple ball is at {target}, inside the wardrobe.")
    print("  - The wardrobe opening is at the front (x=9.5), between y=4 and y=8.")
    # Check from right viewpoint
    v = V_RIGHT
    line_dir = target - v
    # Equation: y = v.y + (target.y - v.y) * (x - v.x) / (target.x - v.x)
    # Find y when x = 9.5 (wardrobe front)
    y_at_opening = v.y + (target.y - v.y) * (9.5 - v.x) / (target.x - v.x)
    print(f"  - From Joe's rightmost viewpoint {v}, the line of sight to the ball intersects")
    print(f"    the wardrobe's front plane (x=9.5) at y = {y_at_opening:.2f}.")
    print("    Equation: y = -0.25 + (4.25 - (-0.25)) * (9.5 - 7.5) / (11.75 - 7.5)")
    if 4 <= y_at_opening <= 8:
        print("  - This is within the opening. The ball might be visible.")
    else:
        print(f"  - Since y={y_at_opening:.2f} is not in the opening range [4, 8], the view is blocked by the wardrobe's front panel.")
    print("  - Visibility: Not Visible\n")
    
    # --- Blue Ball (NE corner) ---
    print("Blue ball (NE corner):")
    target = BALLS['Blue']
    v = V_LEFT # Most advantageous viewpoint
    print(f"  - The Blue ball is at {target}.")
    print(f"  - The main obstacle is the wardrobe's open north door, a plane at y=8 for x in [9.5, 11.5].")
    print(f"  - From the most advantageous viewpoint {v}, we check if the line of sight passes this door.")
    # Find intersection of sightline with the plane y=8
    # Equation: x = v.x + (target.x - v.x) * (y - v.y) / (target.y - v.y)
    line_dir = target - v
    x_at_door_plane = v.x + (target.x - v.x) * (8.0 - v.y) / (target.y - v.y)
    print("  - Equation: x = 4.5 + (11.75 - 4.5) * (8.0 - (-0.25)) / (11.75 - (-0.25))")
    print(f"  - The line of sight intersects the door's plane (y=8) at x = {x_at_door_plane:.3f}.")
    if 9.5 <= x_at_door_plane <= 11.5:
        print(f"  - Since x={x_at_door_plane:.3f} is in the door's range [9.5, 11.5], the view is blocked.")
        print("  - Visibility: Not Visible\n")
    else:
        print(f"  - Since x={x_at_door_plane:.3f} is just in front of the door's range [9.5, 11.5], the sightline just skims past.")
        print("  - Visibility: Visible\n")
        visible_balls.append("Blue")

    # --- Yellow Ball (NW corner) ---
    print("Yellow ball (NW corner):")
    target = BALLS['Yellow']
    v = V_RIGHT # Most advantageous viewpoint
    print(f"  - The Yellow ball is at {target}.")
    print(f"  - The main obstacle is the open room door (a plane at x=4.5 for y in [0,3]).")
    print(f"  - From the widest viewpoint {v}, we check the intersection with the door's plane.")
    # Find intersection of sightline with plane x=4.5
    # Equation: y = v.y + (target.y - v.y) * (x - v.x) / (target.x - v.x)
    line_dir = target - v
    y_at_door_plane = v.y + (target.y - v.y) * (4.5 - v.x) / (target.x - v.x)
    print("  - Equation: y = -0.25 + (11.75 - (-0.25)) * (4.5 - 7.5) / (0.25 - 7.5)")
    print(f"  - The line of sight intersects the door's plane (x=4.5) at y = {y_at_door_plane:.2f}.")
    print(f"  - The physical door occupies y in [0, 3]. Since y={y_at_door_plane:.2f} is above the door, the view is clear.")
    print("  - Visibility: Visible\n")
    visible_balls.append("Yellow")

    # --- Green Ball (on SW shelf) ---
    print("Green ball (on SW shelf):")
    target = BALLS['Green']
    v = V_RIGHT # Most advantageous viewpoint
    print(f"  - The Green ball is at {target}, high up on the shelf.")
    print(f"  - The view is up and to the left. We must check if it's blocked by the top of the doorway (lintel).")
    # Find height of sightline as it passes through the doorway (y=0)
    # Equation: z = v.z + (target.z - v.z) * (y - v.y) / (target.y - v.y)
    line_dir = target - v
    z_at_doorway = v.z + (target.z - v.z) * (0.0 - v.y) / (target.y - v.y)
    print("  - Equation: z = 4.75 + (7.25 - 4.75) * (0.0 - (-0.25)) / (0.25 - (-0.25))")
    print(f"  - The line of sight passes through the doorway opening (y=0) at a height of z = {z_at_doorway:.2f} ft.")
    if z_at_doorway >= 6.67: # Standard door height is 80 inches (6.67 ft)
        print(f"  - This height is above the standard doorway. The view is blocked by the lintel.")
        print("  - Visibility: Not Visible\n")
    else:
        print(f"  - This height is below the standard doorway height of 6.67 ft. The view is clear.")
        print("  - Visibility: Visible\n")
        visible_balls.append("Green")

    print("--- Conclusion ---")
    print("The balls Joe can see are:")
    for ball in sorted(visible_balls):
        print(f"- {ball}")

    # Final answer format
    final_answer = ", ".join(sorted(visible_balls))
    print(f"\n<<<Joe can see the {final_answer} balls.>>>")

solve()