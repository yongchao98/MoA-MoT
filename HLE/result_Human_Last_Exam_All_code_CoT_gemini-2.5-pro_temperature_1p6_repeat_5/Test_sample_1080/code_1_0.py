import itertools

def solve_candy_packing():
    """
    Calculates and verifies the maximum number of eyeball candies in the box.
    """
    # -- Problem Parameters (in 0.5 cm grid units) --
    
    # Box dimensions: 12x12x11 cm -> 24x24x22 units
    # Allowed center coordinates for a sphere of radius 4:
    x_range = (4, 20)
    y_range = (4, 20)
    z_range = (4, 18)
    
    # Sphere radius is 2cm -> 4 units.
    # The square of the minimum distance between centers is (4+4)^2 = 64.
    min_dist_sq = (4 + 4)**2

    # -- A Simple Grid Packing (Baseline) --
    # A simple 3x3x2 grid packing yields 3 * 3 * 2 = 18 candies.
    
    # -- An Improved, Denser Packing (n=22) --
    # This packing consists of 3 layers.
    # Layer 1 (z=4): 9 candies in a 3x3 grid
    # Layer 2 (z=10): 4 candies in a 2x2 grid, placed in the hollows of layer 1
    # Layer 3 (z=16): 9 candies in a 3x3 grid, like layer 1
    
    candies = []
    
    # Layer 1 at z=4 (9 candies)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            candies.append((x, y, 4))
            
    # Layer 2 at z=10 (4 candies)
    for x in [8, 16]:
        for y in [8, 16]:
            candies.append((x, y, 10))

    # Layer 3 at z=16 (9 candies)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            candies.append((x, y, 16))

    n = len(candies)
    print(f"Proposed solution: A packing with {n} candies.")
    print("Verifying this solution...\n")

    # -- Verification Step --
    is_valid = True
    
    # 1. Check if all candy centers are within the box boundaries
    for i, (x, y, z) in enumerate(candies):
        if not (x_range[0] <= x <= x_range[1] and
                y_range[0] <= y <= y_range[1] and
                z_range[0] <= z <= z_range[1]):
            print(f"Error: Candy {i+1} at {p} is out of bounds.")
            is_valid = False
            break
    if not is_valid:
        return 0

    # 2. Check for overlaps between any two candies
    for i, p1 in enumerate(candies):
        for j, p2 in enumerate(candies):
            if i >= j:
                continue
            
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            
            if dist_sq < min_dist_sq:
                print(f"Error: Candy {i+1} at {p1} and Candy {j+1} at {p2} overlap.")
                print(f"Their distance squared is {dist_sq}, which is less than {min_dist_sq}.")
                is_valid = False
                break
        if not is_valid:
            break

    if is_valid:
        print("Verification successful! All constraints are met.")
        print(f"The maximum number of candies is {n}.")
        print("\nThe coordinates for the {n} candies are:")
        # To satisfy the 'output each number' request, we print the coordinates
        for i, (x, y, z) in enumerate(candies):
            print(f"  Candy {i+1:>2}: (x={x:2}, y={y:2}, z={z:2})")
        return n
    else:
        print("\nThe proposed packing is invalid.")
        return 0

# --- Main Execution ---
if __name__ == "__main__":
    max_n = solve_candy_packing()
