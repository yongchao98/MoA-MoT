import math

class Point:
    """A simple class to represent a 3D point."""
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z
    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"

def final_analysis():
    """
    This function performs a step-by-step analysis of which balls Joe can see,
    printing the reasoning for each one.
    """
    # Define key locations and objects
    # Joe is in the middle 3ft of the 12ft south wall, so x is in [4.5, 7.5]
    # He leans in 3 inches (0.25 ft)
    # His eye level is assumed to be 4.75 ft
    joe_vantage_east = Point(7.5, 0.25, 4.75) # Easternmost point in doorway
    joe_vantage_west = Point(4.5, 0.25, 4.75) # Westernmost point in doorway

    # Ball radius is 3 inches (0.25 ft)
    balls = {
        "Red": Point(11.75, 0.25, 0.25),
        "Blue": Point(11.75, 11.75, 2.75), # On a 2.5ft table
        "Yellow": Point(0.25, 11.75, 0.25),
        "Green": Point(0.25, 0.25, 7.25), # On 7ft bookshelf
        "Purple": Point(11.75, 4.25, 0.25)
    }
    
    visible_balls = []

    print("Analyzing Joe's line of sight to each ball:")
    print("="*40)

    # --- 1. Red Ball Analysis ---
    p_joe = joe_vantage_east
    p_ball = balls["Red"]
    print(f"1. Red Ball (SE Corner) at {p_ball}")
    # Line of sight is at y=0.25. Wardrobe starts at y=4.0. No horizontal obstruction.
    m_xz = (p_ball.z - p_joe.z) / (p_ball.x - p_joe.x)
    c_xz = p_joe.z - m_xz * p_joe.x
    print(f"   - View from Joe at {p_joe}.")
    print("   - The line of sight is low along the south wall (y=0.25).")
    print("   - The main potential obstruction, the wardrobe, is located at y=[4, 8].")
    print("   - The view is unobstructed.")
    print(f"   - Equation (x-z plane): z = {m_xz:.3f} * x + {c_xz:.3f}")
    print("   - Verdict: VISIBLE\n")
    visible_balls.append("Red")

    # --- 2. Blue Ball Analysis ---
    p_joe = joe_vantage_east
    p_ball = balls["Blue"]
    print(f"2. Blue Ball (NE Corner) at {p_ball}")
    print(f"   - View from Joe at {p_joe}.")
    # Intersection with wardrobe face at x=9.5
    t = (9.5 - p_joe.x) / (p_ball.x - p_joe.x)
    y_intersect = p_joe.y + t * (p_ball.y - p_joe.y)
    z_intersect = p_joe.z + t * (p_ball.z - p_joe.z)
    print(f"   - The line of sight must pass the wardrobe's front face at x=9.5.")
    print(f"   - The line intersects this plane at y={y_intersect:.2f} and z={z_intersect:.2f}.")
    print("   - The wardrobe body occupies y=[4, 8] and z=[0, 6] (assumed height) at this plane.")
    print(f"   - Since y={y_intersect:.2f} is between 4 and 8, and z={z_intersect:.2f} is between 0 and 6, the view is blocked.")
    print("   - Equation (x-y plane): y - 0.25 = ((11.75-0.25)/(11.75-7.5)) * (x-7.5)")
    print(f"   - y - 0.25 = {(11.5/4.25):.3f} * (x - 7.5)")
    print("   - Verdict: NOT VISIBLE\n")
    
    # --- 3. Yellow Ball Analysis ---
    p_joe = joe_vantage_east
    p_ball = balls["Yellow"]
    print(f"3. Yellow Ball (NW Corner) at {p_ball}")
    print(f"   - View from Joe at {p_joe}.")
    # Intersection with door at x=4.5
    t_door = (4.5 - p_joe.x) / (p_ball.x - p_joe.x)
    y_at_door = p_joe.y + t_door * (p_ball.y - p_joe.y)
    print(f"   - Checking for obstruction by the door (at x=4.5) and bookshelf (at x=1.0).")
    print(f"   - The line of sight intersects the door's plane (x=4.5) at y={y_at_door:.2f}.")
    print("   - This is well outside the door's y-range of [0, 3].")
    print("   - The bookshelf ends at y=4, and the line of sight passes far from it.")
    print(f"   - Equation (x-y plane): y - 0.25 = ((11.75-0.25)/(0.25-7.5)) * (x-7.5)")
    print(f"   - y - 0.25 = {(11.5/-7.25):.3f} * (x - 7.5)")
    print("   - Verdict: VISIBLE\n")
    visible_balls.append("Yellow")

    # --- 4. Green Ball Analysis ---
    p_joe = joe_vantage_east
    p_ball = balls["Green"]
    print(f"4. Green Ball (on SW Bookshelf) at {p_ball}")
    print(f"   - View from Joe at {p_joe}.")
    # Intersection with door at x=4.5
    t = (4.5 - p_joe.x) / (p_ball.x - p_joe.x)
    z_at_door = p_joe.z + t * (p_ball.z - p_joe.z)
    print("   - The view is along the y=0.25 line, so we must check for vertical obstructions.")
    print("   - The main obstruction is the inward-swinging door at x=4.5.")
    print(f"   - The line of sight intersects the door's plane (x=4.5) at height z={z_at_door:.2f} ft.")
    print("   - A standard door is about 6.67 ft tall, so the line of sight is blocked by the door panel.")
    print(f"   - Equation (x-z plane): z - 4.75 = ((7.25-4.75)/(0.25-7.5)) * (x-7.5)")
    print(f"   - z - 4.75 = {(2.5/-7.25):.3f} * (x-7.5)")
    print("   - Verdict: NOT VISIBLE\n")

    # --- 5. Purple Ball Analysis ---
    p_joe = joe_vantage_west # Use west vantage point for best angle
    p_ball = balls["Purple"]
    print(f"5. Purple Ball (in Wardrobe) at {p_ball}")
    print(f"   - View from Joe at {p_joe}.")
    # Intersection with wardrobe opening at x=9.5
    t = (9.5 - p_joe.x) / (p_ball.x - p_joe.x)
    y_intersect = p_joe.y + t * (p_ball.y - p_joe.y)
    print("   - To be seen, the line of sight must pass through the wardrobe's opening at x=9.5, y=[4, 8].")
    print(f"   - From Joe's best viewpoint, the line intersects the x=9.5 plane at y={y_intersect:.2f}.")
    print("   - Since this is lower than the opening's bottom edge (y=4), the view is blocked by the wardrobe's base.")
    print(f"   - Equation to find y-intersect: y_intersect = 0.25 + (({9.5-4.5})/({11.75-4.5})) * (4.25-0.25)")
    print(f"   - y_intersect = 0.25 + (5.0 / 7.25) * 4.0 = {y_intersect:.2f}")
    print("   - Verdict: NOT VISIBLE\n")
    
    print("="*40)
    print("Conclusion:")
    if visible_balls:
        print(f"Joe can see the following balls: {', '.join(visible_balls)}.")
    else:
        print("Joe cannot see any of the balls.")

final_analysis()
<<<Red, Yellow>>>