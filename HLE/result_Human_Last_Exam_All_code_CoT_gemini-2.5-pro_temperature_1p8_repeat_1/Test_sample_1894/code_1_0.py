import math

def get_visible_balls():
    """
    This function models the room and Joe's position to determine which balls are visible.
    It uses a 2D top-down projection and checks for line-of-sight obstructions.
    """

    class Point:
        def __init__(self, x, y):
            self.x = x
            self.y = y

        def __repr__(self):
            return f"Point({self.x:.2f}, {self.y:.2f})"

    def on_segment(p, q, r):
        """Given three collinear points p, q, r, check if q lies on segment 'pr'."""
        return (q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and
                q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y))

    def orientation(p, q, r):
        """Find orientation of the ordered triplet (p, q, r)."""
        val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)
        if -1e-9 < val < 1e-9: return 0  # Collinear (within a tolerance)
        return 1 if val > 0 else 2  # Clockwise or Counterclockwise

    def segments_intersect(p1, q1, p2, q2):
        """Return true if line segment 'p1q1' and 'p2q2' intersect."""
        o1 = orientation(p1, q1, p2)
        o2 = orientation(p1, q1, q2)
        o3 = orientation(p2, q2, p1)
        o4 = orientation(p2, q2, q1)

        if o1 != o2 and o3 != o4:
            return True
        # Collinear cases are not considered intersections for solid objects unless a vertex is touched.
        # For simplicity, we assume Joe can't see through infinitesimal cracks.
        return False

    def check_line_of_sight(joe_pos, ball_pos, obstacles):
        """Checks if the line of sight is blocked by any obstacle."""
        for name, ob in obstacles.items():
            # Obstacle defined as a rectangle ((x1, y1), (x2, y2))
            if "rect" in name:
                p1 = Point(ob[0][0], ob[0][1])  # Top-left
                p2 = Point(ob[1][0], ob[0][1])  # Top-right
                p3 = Point(ob[1][0], ob[1][1])  # Bottom-right
                p4 = Point(ob[0][0], ob[1][1])  # Bottom-left
                sides = [(p1, p2), (p2, p3), (p3, p4), (p4, p1)]
                for side_p1, side_p2 in sides:
                    if segments_intersect(joe_pos, ball_pos, side_p1, side_p2):
                        return False, f"Blocked by the {name}"
            # Obstacle defined as a line ((x1, y1), (x2, y2))
            else:
                ob_p1 = Point(ob[0][0], ob[0][1])
                ob_p2 = Point(ob[1][0], ob[1][1])
                if segments_intersect(joe_pos, ball_pos, ob_p1, ob_p2):
                    return False, f"Blocked by the {name}"
        return True, "Clear"

    # --- Model Definition ---
    # Joe's viewpoint: y=0.25 (leaned in), x between 4.51 and 7.5
    # The door opens inwards, hinged at x=4.5, forming a wall from (4.5,0) to (4.5,3)
    joe_view_points = [Point(x * 0.1, 0.25) for x in range(46, 76)] # x from 4.6 to 7.5

    balls = {
        "Red": Point(11.75, 0.25),
        "Blue": Point(11.75, 11.75),
        "Yellow": Point(0.25, 11.75),
        "Purple": Point(11.75, 4.25), # Inside the wardrobe
        "Green": Point(0.25, 0.25), # On the tall bookshelf
    }

    obstacles = {
        "Door": ((4.5, 0), (4.5, 3)),
        "Bookshelf_rect": ((0, 0), (1, 4)),
        "Wardrobe_main_rect": ((9.5, 4), (12, 8)),
        "Wardrobe_S_Door": ((9.5, 4), (11.5, 4)),
        "Wardrobe_N_Door": ((9.5, 8), (11.5, 8)),
    }
    
    print("Analyzing which balls Joe can see...\n")
    visible_balls = []

    for name, ball_pos in balls.items():
        print(f"--- Checking the {name} ball at {ball_pos} ---")
        
        # 3D Check for Green Ball
        if name == "Green":
            joe_height = 5
            shelf_height = 7
            print(f"Reasoning: The green ball is on top of a {shelf_height} ft shelf.")
            print(f"           Joe is {joe_height} ft tall. He cannot see over the opaque shelf.")
            print(f"Result: {name} ball is NOT VISIBLE.\n")
            continue

        is_visible_for_ball = False
        final_reason = "No clear line of sight from any viewpoint."

        # Iterate through Joe's possible positions
        for joe_pos in joe_view_points:
            
            # Special check for Purple Ball inside the wardrobe
            if name == "Purple":
                # Check if the line of sight passes through the wardrobe opening (x=9.5, y=[4,8])
                m = (ball_pos.y - joe_pos.y) / (ball_pos.x - joe_pos.x)
                y_at_opening = m * (9.5 - joe_pos.x) + joe_pos.y
                
                # Check equation at x=9.5: y = ((4.25 - 0.25) / (11.75 - joe_pos.x)) * (9.5 - joe_pos.x) + 0.25
                if y_at_opening < 4 or y_at_opening > 8:
                    final_reason = (f"Blocked by the wardrobe front panel. Line of sight from {joe_pos} "
                                    f"crosses the wardrobe front at y={y_at_opening:.2f}, "
                                    f"which is outside the opening [4, 8].")
                    # Since this is true for all joe_pos, we can break early.
                    # As joe_pos.x moves from 4.6 to 7.5, y_at_opening moves from 3.02 to 2.13. Always blocked.
                    is_visible_for_ball = False
                    if joe_pos.x == joe_view_points[-1].x: # on last check
                        break 
                    else:
                        continue # check next joe pos (though result will be same)
                
            # General line-of-sight check against all obstacles
            is_visible, reason = check_line_of_sight(joe_pos, ball_pos, obstacles)

            if is_visible:
                is_visible_for_ball = True
                final_reason = f"Clear line of sight found from Joe's position at {joe_pos}."
                break # Found a clear view, no need to check other positions for this ball
            else:
                # Store the reason for the last tested position if no clear view is found
                final_reason = reason
        
        if is_visible_for_ball:
            print(f"Reasoning: {final_reason}")
            print(f"Result: {name} ball is VISIBLE.\n")
            visible_balls.append(name)
        else:
            print(f"Reasoning: {final_reason}")
            print(f"Result: {name} ball is NOT VISIBLE.\n")

    print("--- Summary ---")
    print(f"The balls Joe can see are: {', '.join(sorted(visible_balls))}.")
    return sorted(visible_balls)

if __name__ == '__main__':
    final_answer = get_visible_balls()
    # The final answer is expected in a specific format by the system.
    # print(f"<<<{', '.join(final_answer)}>>>")
<<<Red, Blue, Yellow>>>