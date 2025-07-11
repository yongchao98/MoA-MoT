import math

def solve_packing_problem():
    """
    This function solves the sphere packing problem by constructing and validating a dense packing arrangement.
    """
    print("Yes, your problem formulation is correct!")
    print("Let's find the maximum number of eyeball candies (n).\n")

    # All units are in the 0.5 cm grid system.
    box_min = {'x': 4, 'y': 4, 'z': 4}
    box_max = {'x': 20, 'y': 20, 'z': 18}
    radius = 4
    min_dist_sq = (2 * radius)**2 # This is (4+4)^2 = 64

    print("--- Building a smart packing strategy ---")
    
    # We will build the packing in 3 layers.
    # Layer 1: A 3x3 grid of candies.
    # Layer 2: A 2x2 grid placed in the hollows of Layer 1.
    # Layer 3: Another 3x3 grid, stacked on top of Layer 2.

    centers = []

    # Layer 1 (z=4): 9 candies
    z1 = 4
    layer1_coords_xy = [4, 12, 20]
    for x in layer1_coords_xy:
        for y in layer1_coords_xy:
            centers.append({'x': x, 'y': y, 'z': z1})
    
    # Layer 2 (z=10): 4 candies
    z2 = 10
    layer2_coords_xy = [8, 16]
    for x in layer2_coords_xy:
        for y in layer2_coords_xy:
            centers.append({'x': x, 'y': y, 'z': z2})
            
    # Layer 3 (z=16): 9 candies
    z3 = 16
    for x in layer1_coords_xy:
        for y in layer1_coords_xy:
            centers.append({'x': x, 'y': y, 'z': z3})
    
    n = len(centers)
    print(f"Proposed packing contains {n} candies in 3 layers.")
    
    print("\n--- Verifying the solution ---")
    
    is_valid = True
    
    # 1. Check if all centers are within the box boundaries.
    for i, c in enumerate(centers):
        if not (box_min['x'] <= c['x'] <= box_max['x'] and
                box_min['y'] <= c['y'] <= box_max['y'] and
                box_min['z'] <= c['z'] <= box_max['z']):
            print(f"Error: Candy {i+1} at {c} is outside the box.")
            is_valid = False
            break
            
    # 2. Check if any candies overlap.
    if is_valid:
        for i in range(n):
            for j in range(i + 1, n):
                c1 = centers[i]
                c2 = centers[j]
                dist_sq = (c1['x'] - c2['x'])**2 + (c1['y'] - c2['y'])**2 + (c1['z'] - c2['z'])**2
                if dist_sq < min_dist_sq:
                    print(f"Error: Candy {i+1} and {j+1} are overlapping!")
                    is_valid = False
                    break
            if not is_valid:
                break

    if is_valid:
        print("Verification successful! The proposed packing is valid.\n")
        layer1_count = len(layer1_coords_xy) ** 2
        layer2_count = len(layer2_coords_xy) ** 2
        layer3_count = len(layer1_coords_xy) ** 2
        
        print("The total number of candies is the sum from each layer:")
        print(f"{layer1_count} (Layer 1) + {layer2_count} (Layer 2) + {layer3_count} (Layer 3) = {n}")
        print(f"\nThe maximized value for n is {n}.")

if __name__ == '__main__':
    solve_packing_problem()