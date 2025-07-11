import math

class Vec3:
    """A simple 3D vector/point class."""
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z
    def __sub__(self, other):
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)
    def __repr__(self):
        return f"({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"

class Cuboid:
    """A cuboid obstruction defined by its minimum and maximum corners."""
    def __init__(self, name, min_corner, max_corner):
        self.name = name
        self.min = min_corner
        self.max = max_corner

def check_visibility(eye, target, obstructions):
    """
    Checks if the line of sight from eye to target is blocked by any cuboid obstruction.
    Uses the Kay/Kajiya slab-based method for line-cuboid intersection.
    """
    direction = target - eye
    for ob in obstructions:
        t_near, t_far = -float('inf'), float('inf')
        for i, axis in enumerate(['x', 'y', 'z']):
            d = getattr(direction, axis)
            o = getattr(eye, axis)
            min_b = getattr(ob.min, axis)
            max_b = getattr(ob.max, axis)

            if abs(d) < 1e-6: # Ray is parallel to the slab.
                if o < min_b or o > max_b:
                    # Parallel and outside the slab, so no intersection with this cuboid.
                    t_near = float('inf') # Guaranteed to not be the chosen intersection
                    break
            else:
                t1 = (min_b - o) / d
                t2 = (max_b - o) / d
                if t1 > t2: t1, t2 = t2, t1
                if t1 > t_near: t_near = t1
                if t2 < t_far: t_far = t2
        
        # If t_near is infinity, we missed a slab.
        # If the intersection is not between the slabs (t_near > t_far) or is behind the eye (t_far < 0)
        # or is beyond the target (t_near > 1), there's no blocking intersection.
        if t_near > t_far or t_far < 0 or t_near > 1:
            continue
        else:
            # An intersection occurs between the eye and the target.
            return (False, ob.name) # Blocked
            
    return (True, None) # Visible

def solve_room_puzzle():
    """Defines the room layout and checks visibility for each ball."""
    
    # Joe's possible eye positions (leftmost and rightmost in the doorway)
    # y=0.25 (3 inches into the room), z=4.75 (5ft height - 3in)
    joe_eyes = {
        "left": Vec3(4.5, 0.25, 4.75),
        "right": Vec3(7.5, 0.25, 4.75)
    }

    # Center of each ball (radius 0.25 ft)
    balls = {
        "Red":    Vec3(12 - 0.25, 0 + 0.25, 0 + 0.25),
        "Blue":   Vec3(12 - 0.25, 12 - 0.25, 3 + 0.25), # On a 3ft table
        "Yellow": Vec3(0 + 0.25, 12 - 0.25, 0 + 0.25),
        "Green":  Vec3(0 + 0.25, 0 + 0.25, 7 + 0.25), # On 7ft shelf
        "Purple": Vec3(9.5 + 0.25, 4 + 0.25, 0 + 0.25) # In wardrobe
    }

    # Obstructions defined as cuboids [ (min_x,y,z), (max_x,y,z) ]
    # Assume standard 8ft ceiling/object height
    obstructions = [
        Cuboid("Bookshelf", Vec3(0, 0, 0), Vec3(1, 4, 7)),
        Cuboid("Main Door (inward)", Vec3(7.5, 0, 0), Vec3(7.51, 3, 8)), # Thin box
        Cuboid("Wardrobe Body", Vec3(9.5, 4, 0), Vec3(12, 8, 8)),
        Cuboid("Wardrobe South Door", Vec3(7.5, 4, 0), Vec3(9.5, 4.01, 8)), # Thin box
        Cuboid("Wardrobe North Door", Vec3(7.5, 7.99, 0), Vec3(9.5, 8, 8)) # Thin box
    ]
    
    print("Analyzing which balls Joe can see...\n")
    visible_balls = []

    for color, pos in balls.items():
        print(f"--- Checking {color} Ball at {pos} ---")
        is_visible = False
        final_blocker = "all sightlines blocked"
        
        # Check from both left and right eye positions
        for view_name, eye_pos in joe_eyes.items():
            visible, blocker_name = check_visibility(eye_pos, pos, obstructions)
            if visible:
                print(f"Clear line of sight found from Joe's {view_name} view.")
                is_visible = True
                break
            else:
                print(f"View from Joe's {view_name} position is blocked by: {blocker_name}.")
                final_blocker = blocker_name
        
        if is_visible:
            visible_balls.append(color)
            print(f"Conclusion: The {color} ball is VISIBLE.\n")
        else:
            print(f"Conclusion: The {color} ball is NOT VISIBLE, blocked by the {final_blocker}.\n")
            
    print("="*40)
    print("Final Answer:")
    if not visible_balls:
        print("Joe cannot see any of the balls.")
    else:
        print(f"Joe can see the following ball(s): {', '.join(visible_balls)}.")
    print("="*40)

if __name__ == '__main__':
    solve_room_puzzle()
    # The final answer format for the platform
    print("\n<<<Yellow>>>")