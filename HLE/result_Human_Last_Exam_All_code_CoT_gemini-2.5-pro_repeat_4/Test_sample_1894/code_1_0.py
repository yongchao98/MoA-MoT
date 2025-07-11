import math

def check_visibility():
    """
    Analyzes the visibility of balls in a room based on a 2D floor plan.
    """

    # --- 1. Define Objects and Balls ---

    # Joe's viewing position is a line segment he can move along
    # He can lean in 3 inches (0.25 feet)
    joe_view_start = (4.5, 0.25)
    joe_view_end = (7.5, 0.25)

    # Ball positions (using the center of the 6-inch diameter balls)
    balls = {
        "Red": (11.75, 0.25),
        "Blue": (11.75, 11.75),
        "Yellow": (0.25, 11.75)
    }

    # Obstructions are defined as line segments
    obstructions = []
    # Open room door (hinged on the east side of the doorway)
    obstructions.append(((7.5, 0), (7.5, 3)))
    # Bookshelf (1ft deep, 4ft wide along west wall)
    bookshelf_rect = [(0, 0), (1, 0), (1, 4), (0, 4)]
    obstructions.append((bookshelf_rect[0], bookshelf_rect[1]))
    obstructions.append((bookshelf_rect[1], bookshelf_rect[2]))
    obstructions.append((bookshelf_rect[2], bookshelf_rect[3]))
    obstructions.append((bookshelf_rect[3], bookshelf_rect[0]))
    # Wardrobe (2.5ft deep, 4ft wide on east wall)
    wardrobe_rect = [(9.5, 4), (12, 4), (12, 8), (9.5, 8)]
    obstructions.append((wardrobe_rect[0], wardrobe_rect[1]))
    obstructions.append((wardrobe_rect[1], wardrobe_rect[2]))
    obstructions.append((wardrobe_rect[2], wardrobe_rect[3]))
    obstructions.append((wardrobe_rect[3], wardrobe_rect[0]))
    # Wardrobe doors (open 90 degrees)
    obstructions.append(((7.5, 4), (9.5, 4))) # South door
    obstructions.append(((7.5, 8), (9.5, 8))) # North door

    # --- 2. Define Intersection Logic ---

    def intersects(los_start, los_end, obs_start, obs_end):
        """
        Checks if the line of sight (los) intersects an obstruction line segment.
        Returns True if they intersect, False otherwise.
        """
        p0_x, p0_y = los_start
        p1_x, p1_y = los_end
        p2_x, p2_y = obs_start
        p3_x, p3_y = obs_end

        s1_x = p1_x - p0_x
        s1_y = p1_y - p0_y
        s2_x = p3_x - p2_x
        s2_y = p3_y - p2_y

        # Determinant
        denominator = (-s2_x * s1_y + s1_x * s2_y)
        if denominator == 0:
            return False # Parallel lines

        # Solve for t and u
        s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / denominator
        t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / denominator

        # An intersection exists if 0 < t < 1 and 0 < s < 1
        # We use a small epsilon to avoid floating point issues at endpoints
        epsilon = 1e-9
        if (epsilon < s < 1 - epsilon) and (epsilon < t < 1 - epsilon):
            return True
        return False

    # --- 3. Perform Visibility Check ---
    
    print("Analyzing visibility from Joe's perspective...")
    print("-" * 30)

    # Initial analysis of special cases
    print("Green ball: NOT VISIBLE. It is on a 7-foot-tall shelf, which is above Joe's eye level (he is 5 feet tall).")
    print("Purple ball: NOT VISIBLE. It is inside the wardrobe, occluded by the wardrobe's own front and side panels.")
    print("-" * 30)

    visible_balls = []

    # Check remaining balls against all obstructions
    for ball_name, ball_pos in balls.items():
        is_visible = False
        # Check from 100 different points along Joe's viewing line for thoroughness
        for i in range(101):
            t = i / 100.0
            joe_x = joe_view_start[0] * (1 - t) + joe_view_end[0] * t
            joe_y = joe_view_start[1] * (1 - t) + joe_view_end[1] * t
            joe_pos = (joe_x, joe_y)

            view_is_blocked = False
            for obs in obstructions:
                if intersects(joe_pos, ball_pos, obs[0], obs[1]):
                    view_is_blocked = True
                    break # This line of sight is blocked, try next position for Joe
            
            if not view_is_blocked:
                is_visible = True
                break # Found a clear line of sight, no need to check other positions
        
        if is_visible:
            print(f"{ball_name} ball: VISIBLE. A clear line of sight was found.")
            visible_balls.append(ball_name)
        else:
            print(f"{ball_name} ball: NOT VISIBLE. The line of sight is blocked by an obstruction.")

    print("-" * 30)
    print("Final Conclusion:")
    if visible_balls:
        print(f"Joe can see the following ball(s): {', '.join(visible_balls)}.")
    else:
        print("Joe cannot see any of the balls.")


check_visibility()
<<<Yellow>>>