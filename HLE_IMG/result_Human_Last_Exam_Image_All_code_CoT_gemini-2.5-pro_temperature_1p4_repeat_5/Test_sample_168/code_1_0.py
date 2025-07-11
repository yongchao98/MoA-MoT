import math

def solve_mirror_paths():
    """
    Calculates the number of light ray paths given the number of reflections on each mirror.
    """
    # Step 1: Define the number of reflections for each mirror from the problem statement.
    # m = 2 reflections on G1.
    n_g1 = 2  # Mirror G1 at x=0 (left)
    n_g2 = 1  # Mirror G2 at y=0 (bottom)
    n_g3 = 1  # Mirror G3 at x=w (right)
    n_g4 = 1  # Mirror G4 at y=h (top)

    print("Step 1: Define reflection counts for each mirror.")
    print(f"Reflections on G1 (left, m=2): {n_g1}")
    print(f"Reflections on G2 (bottom): {n_g2}")
    print(f"Reflections on G3 (right): {n_g3}")
    print(f"Reflections on G4 (top): {n_g4}")
    print("-" * 40)

    # Step 2: Calculate total horizontal and vertical reflections.
    n_h = n_g1 + n_g3
    n_v = n_g2 + n_g4

    print("Step 2: Calculate total horizontal and vertical reflections.")
    print(f"Total horizontal reflections (G1+G3): {n_g1} + {n_g3} = {n_h}")
    print(f"Total vertical reflections (G2+G4): {n_g2} + {n_g4} = {n_v}")
    print("-" * 40)
    
    # Step 3 & 4: Determine valid path directions and corresponding image floor values.
    possible_floor_pairs = []
    print("Step 3: Determine which path directions satisfy the reflection counts.")

    # A path can be left/right and up/down. We check the four combinations.
    
    # Check horizontal direction (left vs. right)
    # Leftward path: G1 reflections = ceil(n_h/2), G3 reflections = floor(n_h/2)
    if math.ceil(n_h / 2) == n_g1 and math.floor(n_h / 2) == n_g3:
        is_leftward = True
        Fx = -n_h
        print(f"- Horizontal path must be LEFTWARD (G1 hits={math.ceil(n_h / 2)}, G3 hits={math.floor(n_h / 2)}).")
    # Rightward path: G1 reflections = floor(n_h/2), G3 reflections = ceil(n_h/2)
    elif math.floor(n_h / 2) == n_g1 and math.ceil(n_h / 2) == n_g3:
        is_leftward = False
        Fx = n_h
        print(f"- Horizontal path must be RIGHTWARD (G1 hits={math.floor(n_h / 2)}, G3 hits={math.ceil(n_h / 2)}).")
    else: # Should not happen with given numbers
        is_leftward = None
        
    # Check vertical direction (up vs. down)
    ways_from_direction = []
    
    # Check upward path
    # G2 reflections = floor(n_v/2), G4 reflections = ceil(n_v/2)
    if math.floor(n_v / 2) == n_g2 and math.ceil(n_v / 2) == n_g4:
        Fy = n_v
        print(f"- A valid vertical path is UPWARD (G2 hits={math.floor(n_v / 2)}, G4 hits={math.ceil(n_v / 2)}).")
        if is_leftward is not None:
             possible_floor_pairs.append({'direction': 'Left-Up', 'Fx': Fx, 'Fy': Fy, 'count': 1})

    # Check downward path
    # G2 reflections = ceil(n_v/2), G4 reflections = floor(n_v/2)
    if math.ceil(n_v / 2) == n_g2 and math.floor(n_v / 2) == n_g4:
        Fy = -n_v
        print(f"- A valid vertical path is DOWNWARD (G2 hits={math.ceil(n_v / 2)}, G4 hits={math.floor(n_v / 2)}).")
        if is_leftward is not None:
            possible_floor_pairs.append({'direction': 'Left-Down', 'Fx': Fx, 'Fy': Fy, 'count': 1})
            
    print("-" * 40)

    # Step 5: Sum the number of ways. Each valid direction corresponds to one family of images, hence one way.
    print("Step 4: Sum the ways from each valid path direction.")
    total_ways = sum(p['count'] for p in possible_floor_pairs)
    
    ways_list = [str(p['count']) for p in possible_floor_pairs]
    equation_str = " + ".join(ways_list)

    for p in possible_floor_pairs:
        print(f"Path direction {p['direction']} provides {p['count']} way.")

    print(f"\nThe total number of ways is {equation_str} = {total_ways}.")
    
    return total_ways

solve_mirror_paths()
<<<2>>>