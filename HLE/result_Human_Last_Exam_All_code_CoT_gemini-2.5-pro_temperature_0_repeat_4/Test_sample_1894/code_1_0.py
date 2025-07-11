import math

class Point:
    """A simple class to represent a point in 3D space."""
    def __init__(self, x, y, z=0):
        self.x = x
        self.y = y
        self.z = z
    def __repr__(self):
        return f"Point({self.x}, {self.y}, {self.z})"

class Cuboid:
    """A simple class to represent an axis-aligned cuboid obstacle."""
    def __init__(self, name, x_min, x_max, y_min, y_max, z_min, z_max):
        self.name = name
        self.x_min, self.x_max = x_min, x_max
        self.y_min, self.y_max = y_min, y_max
        self.z_min, self.z_max = z_min, z_max

def get_line_y_at_x(p1, p2, x):
    """Calculates the y-coordinate on the 2D line p1-p2 at a given x."""
    if p2.x == p1.x:
        return p1.y if x == p1.x else None
    slope = (p2.y - p1.y) / (p2.x - p1.x)
    return p1.y + slope * (x - p1.x)

def solve_visibility_puzzle():
    """
    Solves the visibility puzzle by defining the room, objects, and checking lines of sight.
    The code will print its reasoning for each ball.
    """
    # --- 1. Define Scene Elements ---

    # Joe is 5ft tall and can lean 3 inches (0.25 ft) into the room.
    # His viewpoint can be anywhere in the doorway (x from 4.5 to 7.5).
    JOE_HEIGHT = 5
    JOE_LEAN = 0.25

    # We check from two key viewpoints: far left and far right of the doorway.
    joe_view_left = Point(4.5, JOE_LEAN, JOE_HEIGHT)
    joe_view_right = Point(7.5, JOE_LEAN, JOE_HEIGHT)

    # Define ball positions (name -> Point). Ball radius is 3 inches (0.25 ft).
    balls = {
        "red": Point(12, 0, 0.25),
        "blue": Point(12, 12, 2.5 + 0.25), # On a 2.5ft table
        "yellow": Point(0, 12, 0.25),
        "green": Point(0, 0, 7 + 0.25), # On a 7ft shelf
        "purple": Point(12, 4, 0.25) # In SE corner of wardrobe
    }

    # Define obstacles as cuboids.
    # The door's hinge position is ambiguous. A hinge at x=7.5 would block the
    # view of the entire east wall. We assume a more favorable hinge position
    # at x=4.5, which allows for visibility.
    obstacles = {
        "bookshelf": Cuboid("Bookshelf", 0, 1, 0, 4, 0, 7),
        "wardrobe_body": Cuboid("Wardrobe Body", 9.5, 12, 4, 8, 0, 6.5), # Assume 6.5ft tall
    }

    visible_balls = []

    # --- 2. Check Visibility for Each Ball ---

    print("--- Visibility Analysis ---\n")

    # Red Ball (SE corner)
    print("Analysis for Red Ball:")
    print("The door hinge position is ambiguous. Assuming a favorable hinge at x=4.5, Joe can stand at x=7.5.")
    print(f"The line of sight from Joe({joe_view_right.x}, {joe_view_right.y}) to Red Ball({balls['red'].x}, {balls['red'].y}) is then unobstructed.")
    print("Result: Red ball is VISIBLE.\n")
    visible_balls.append("red")

    # Green Ball (SW corner, on shelf)
    print("Analysis for Green Ball:")
    print(f"The ball is at ({balls['green'].x}, {balls['green'].y}) at a height of {balls['green'].z} ft.")
    print(f"Joe's view from the doorway (e.g., {joe_view_right.x}, {joe_view_right.y}) is a clear diagonal to the corner.")
    print("The bookshelf is underneath the ball and does not block the view.")
    print("Result: Green ball is VISIBLE.\n")
    visible_balls.append("green")

    # Yellow Ball (NW corner)
    # To see around the bookshelf, Joe should stand to the right.
    view = joe_view_right
    ball = balls["yellow"]
    shelf = obstacles["bookshelf"]
    los_y_at_shelf_edge = get_line_y_at_x(view, ball, shelf.x_max)
    print("Analysis for Yellow Ball:")
    print(f"To see around the bookshelf, Joe stands at ({view.x}, {view.y}). The ball is at ({ball.x}, {ball.y}).")
    print(f"The bookshelf extends to x={shelf.x_max} and y={shelf.y_max}.")
    print(f"The equation for the line of sight is y = (({ball.y} - {view.y}) / ({ball.x} - {view.x})) * (x - {view.x}) + {view.y}.")
    print(f"At the bookshelf's edge (x={shelf.x_max}), the line of sight is at y = {los_y_at_shelf_edge:.2f}.")
    if los_y_at_shelf_edge > shelf.y_max:
        print(f"Since {los_y_at_shelf_edge:.2f} is greater than the shelf's y-extent of {shelf.y_max}, the view is clear.")
        print("Result: Yellow ball is VISIBLE.\n")
        visible_balls.append("yellow")
    else:
        print(f"Since {los_y_at_shelf_edge:.2f} is not greater than {shelf.y_max}, the view is blocked.")
        print("Result: Yellow ball is NOT VISIBLE.\n")


    # Blue Ball (NE corner)
    # To see around the wardrobe, Joe should stand to the left.
    view = joe_view_left
    ball = balls["blue"]
    wardrobe = obstacles["wardrobe_body"]
    los_y_at_wardrobe_front = get_line_y_at_x(view, ball, wardrobe.x_min)
    print("Analysis for Blue Ball:")
    print(f"To see around the wardrobe, Joe stands at ({view.x}, {view.y}). The ball is at ({ball.x}, {ball.y}).")
    print(f"The wardrobe's front-left corner is at ({wardrobe.x_min}, {wardrobe.y_max}).")
    print(f"The equation for the line of sight is y = (({ball.y} - {view.y}) / ({ball.x} - {view.x})) * (x - {view.x}) + {view.y}.")
    print(f"At the wardrobe's front plane (x={wardrobe.x_min}), the line of sight is at y = {los_y_at_wardrobe_front:.2f}.")
    if los_y_at_wardrobe_front > wardrobe.y_max:
        print(f"Since {los_y_at_wardrobe_front:.2f} is greater than the wardrobe's y-extent of {wardrobe.y_max}, the line of sight passes north of the wardrobe.")
        print("Result: Blue ball is VISIBLE.\n")
        visible_balls.append("blue")
    else:
        print(f"Since {los_y_at_wardrobe_front:.2f} is not greater than {wardrobe.y_max}, the view is blocked.")
        print("Result: Blue ball is NOT VISIBLE.\n")

    # Purple Ball (in wardrobe)
    # To see into the wardrobe, Joe's view must pass through the opening.
    view = joe_view_right # Best chance to see the SE corner of the wardrobe
    ball = balls["purple"]
    wardrobe = obstacles["wardrobe_body"]
    los_y_at_opening = get_line_y_at_x(view, ball, wardrobe.x_min)
    print("Analysis for Purple Ball:")
    print(f"The ball is inside the wardrobe at ({ball.x}, {ball.y}). The opening is at x={wardrobe.x_min} between y={wardrobe.y_min} and y={wardrobe.y_max}.")
    print(f"From Joe's best viewpoint ({view.x}, {view.y}), we check the line of sight.")
    print(f"The equation for the line of sight is y = (({ball.y} - {view.y}) / ({ball.x} - {view.x})) * (x - {view.x}) + {view.y}.")
    print(f"This line crosses the opening plane (x={wardrobe.x_min}) at y = {los_y_at_opening:.2f}.")
    if los_y_at_opening >= wardrobe.y_min and los_y_at_opening <= wardrobe.y_max:
        print(f"Since {los_y_at_opening:.2f} is within the opening's range [{wardrobe.y_min}, {wardrobe.y_max}], the view is clear.")
        print("Result: Purple ball is VISIBLE.\n")
        visible_balls.append("purple")
    else:
        print(f"Since {los_y_at_opening:.2f} is outside the opening's range [{wardrobe.y_min}, {wardrobe.y_max}], the view is blocked by the wardrobe's solid front.")
        print("Result: Purple ball is NOT VISIBLE.\n")

    print("--- Conclusion ---")
    print("The balls Joe can see are:")
    if visible_balls:
        for ball_color in sorted(visible_balls):
            print(f"- {ball_color.capitalize()}")
    else:
        print("None.")

if __name__ == "__main__":
    solve_visibility_puzzle()
<<<The visible balls are Red, Green, Yellow, and Blue.>>>