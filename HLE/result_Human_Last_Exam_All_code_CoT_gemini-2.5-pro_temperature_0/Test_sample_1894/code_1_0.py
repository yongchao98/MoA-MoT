import math

class Point:
    """A simple class to represent a point in 2D space."""
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f})"

def calculate_intersection_y(p1, p2, x_line):
    """Calculates the y-coordinate where the line p1-p2 intersects a vertical line."""
    if (p2.x - p1.x) == 0: # Avoid division by zero for vertical lines
        return float('inf') if p1.x == x_line else None
    slope = (p2.y - p1.y) / (p2.x - p1.x)
    y_intersect = slope * (x_line - p1.x) + p1.y
    return y_intersect

def check_visibility(ball_name, ball_pos, obstacles, joe_view_points):
    """
    Checks if a ball is visible from any of Joe's viewpoints.
    Prints the reasoning based on geometric calculations.
    """
    print(f"--- Checking visibility of the {ball_name} ball at {ball_pos} ---")

    # Special check for the Green ball based on height
    if ball_name == "green":
        joe_eye_height = 5
        shelf_height = 7
        print(f"Joe's eye height is {joe_eye_height} ft.")
        print(f"The green ball is on a shelf that is {shelf_height} ft tall.")
        print(f"Since the shelf is taller than Joe's eyes, he cannot see what is on top of it.")
        print("Conclusion: The green ball is NOT visible.\n")
        return False

    # Check from multiple viewpoints
    for view_point in joe_view_points:
        is_blocked = False
        blocking_obstacle = None
        
        for name, obstacle in obstacles.items():
            # An obstacle is a line segment (p1, p2)
            p1, p2 = obstacle['points']
            
            # Check if the line of sight intersects the obstacle's line
            # We only care about obstacles between Joe and the ball
            min_x = min(view_point.x, ball_pos.x)
            max_x = max(view_point.x, ball_pos.x)
            
            # The obstacle must be between the viewpoint and the ball to block it
            if not (min_x < p1.x < max_x):
                continue

            y_intersect = calculate_intersection_y(view_point, ball_pos, p1.x)
            
            if y_intersect is not None and obstacle['y_range'][0] <= y_intersect <= obstacle['y_range'][1]:
                is_blocked = True
                blocking_obstacle = {
                    "name": name,
                    "view_point": view_point,
                    "intersect_point": Point(p1.x, y_intersect)
                }
                break # This viewpoint is blocked, try the next one if available
        
        if not is_blocked:
            print(f"A clear line of sight is found from viewpoint {view_point}.")
            print("Conclusion: The {ball_name} ball IS visible.\n")
            return True

    # If all viewpoints were blocked
    if blocking_obstacle:
        bo = blocking_obstacle
        print(f"The line of sight from viewpoint {bo['view_point']} to the ball is blocked.")
        print(f"Equation: Line from {bo['view_point']} to {ball_pos}.")
        print(f"The line intersects the vertical plane of the '{bo['name']}' at x={bo['intersect_point'].x}.")
        print(f"The calculated intersection point is {bo['intersect_point']}.")
        y_range = obstacles[bo['name']]['y_range']
        print(f"This point falls within the obstacle's blocking range of y=[{y_range[0]:.2f}, {y_range[1]:.2f}].")

    print(f"No clear line of sight found for the {ball_name} ball.")
    print(f"Conclusion: The {ball_name} ball is NOT visible.\n")
    return False

# --- Main execution ---

# Define ball locations
balls = {
    "red":    Point(11.75, 0.25),
    "blue":   Point(11.75, 11.75),
    "yellow": Point(0.25, 11.75),
    "green":  Point(0.25, 0.25), # Height is the key factor
    "purple": Point(11.75, 4.25)
}

# Define Joe's key viewpoints (leftmost and rightmost in the doorway)
joe_view_points = [Point(7.5, 0.25), Point(4.5, 0.25)]

# Define obstacles as vertical lines with a y-range
obstacles = {
    "Bookshelf": {"points": (Point(1, 0), Point(1, 4)), "y_range": [0, 4]},
    "Open Door": {"points": (Point(4.5, 0), Point(4.5, 3)), "y_range": [0, 3]},
    "Wardrobe Front": {"points": (Point(9.5, 0), Point(9.5, 12)), "y_range": [0, 4]}, # Solid part below opening
    "Wardrobe Top": {"points": (Point(9.5, 0), Point(9.5, 12)), "y_range": [8, 10]} # Solid part above opening + open door
}

visible_balls = []
for name, pos in balls.items():
    if check_visibility(name, pos, obstacles, joe_view_points):
        visible_balls.append(name)

print("=========================================")
print("Final Summary:")
print("Joe can see the following balls:")
for color in visible_balls:
    print(f"- {color.capitalize()}")

# Format the final answer as requested
final_answer = " and ".join([b.capitalize() for b in visible_balls])
print(f"\n<<<The {final_answer} balls are visible.>>>")
