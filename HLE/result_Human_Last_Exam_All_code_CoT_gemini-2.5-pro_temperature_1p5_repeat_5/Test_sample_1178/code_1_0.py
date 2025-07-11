import math

def get_divisors(n):
    """Find all pairs of divisors (L, W) for a number n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add((i, n//i))
    return list(divs)

def is_formable(n, sides):
    """Check if n can be formed by a non-negative integer sum of side lengths."""
    # Using dynamic programming (unbounded knapsack / change-making problem)
    dp = [False] * (n + 1)
    dp[0] = True
    for i in range(1, n + 1):
        for s in sides:
            if i >= s and dp[i - s]:
                dp[i] = True
                break
    return dp[n]

def solve():
    """
    Find the smallest area of a rectangle admitting a non-guillotine tiling
    with squares from the set S={2x2, 3x3, 5x5, 7x7}.
    """
    SIDES = {2, 3, 5, 7}
    min_area = float('inf')
    best_config = {}

    print("Step 1: The key to a non-guillotine tiling is a 'pinwheel' of four squares.")
    print("This pinwheel is made of two s1 x s1 squares and two s2 x s2 squares.")
    print("This arrangement creates a central hole of size |s1 - s2| x |s1 - s2|.")
    print("To create a complete tiling, this hole must be filled with another square from our set S.\n")

    print("Step 2: Find pairs (s1, s2) from {2, 3, 5, 7} where s_hole = |s1 - s2| is also in the set.")
    
    candidate_areas = []

    # Generate unique pairs from SIDES
    side_list = sorted(list(SIDES))
    for i in range(len(side_list)):
        for j in range(i + 1, len(side_list)):
            s1 = side_list[i]
            s2 = side_list[j]
            s_hole = abs(s1 - s2)

            if s_hole in SIDES:
                area = 2 * s1**2 + 2 * s2**2 + s_hole**2
                tile_sides = {s1, s2, s_hole}
                candidate_areas.append({'s1': s1, 's2': s2, 's_hole': s_hole, 'area': area, 'tile_sides': tile_sides})
    
    print(f"Found {len(candidate_areas)} candidate tiling configurations.\n")

    print("Step 3: For each candidate, calculate the total area and check if a rectangle with this area can be formed.")
    for config in sorted(candidate_areas, key=lambda x: x['area']):
        area = config['area']
        tile_sides = config['tile_sides']
        
        print(f"\nChecking configuration with tiles {config['s1']}x{config['s1']}, {config['s2']}x{config['s2']}, and hole-filler {config['s_hole']}x{config['s_hole']}.")
        print(f"Total area = 2*{config['s1']}^2 + 2*{config['s2']}^2 + {config['s_hole']}^2 = {area}")

        is_viable = False
        divs = get_divisors(area)
        for L, W in divs:
            if is_formable(L, tile_sides) and is_formable(W, tile_sides):
                print(f"Found viable rectangle dimensions: {L} x {W}.")
                is_viable = True
                break
        
        if is_viable:
            if area < min_area:
                min_area = area
                best_config = config
                print(f"This is the smallest viable area found so far.")
        else:
            print(f"No viable rectangle dimensions found for area {area} (e.g., for prime area, a 1x{area} rectangle cannot be tiled).")

    print("\nStep 4: Conclusion")
    print("The smallest valid area corresponds to the first viable configuration found.")
    s1 = best_config['s1']
    s2 = best_config['s2']
    s_hole = best_config['s_hole']
    
    print(f"\nThe smallest integer length rectangle admitting such a tiling is found using squares of side lengths {s1} and {s2} to form the pinwheel, which creates a hole of side length {s_hole}.")
    print("The rectangle is composed of two " + str(s1) + "x" + str(s1) + " squares, two " + str(s2) + "x" + str(s2) + " squares, and one " + str(s_hole) + "x" + str(s_hole) + " square.")
    print("The final area calculation is:")
    print(f"Area = 2 * {s1}^2 + 2 * {s2}^2 + {s_hole}^2 = 2 * {s1*s1} + 2 * {s2*s2} + {s_hole*s_hole} = {best_config['area']}")
    print(f"The area of this rectangle is {min_area}.")

solve()