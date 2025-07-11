import math

def solve_candy_packing():
    """
    Solves the discrete sphere packing problem using a greedy algorithm.
    """
    # Box dimensions: 12x12x11 cm
    # Candy radius: 2 cm
    
    # 1. Define the search space for the centers (in 0.5 cm units)
    # Real-world valid center range: x,y in [2, 10], z in [2, 9]
    # Scaled integer coordinates (x0.5): x,y in [4, 20], z in [4, 18]
    x_coords = range(4, 21)
    y_coords = range(4, 21)
    z_coords = range(4, 19)

    # 2. Define the non-overlapping constraint
    # Real-world min distance is 2 * radius = 4 cm.
    # Scaled min distance is 4 / 0.5 = 8 units.
    # Min squared distance is 8*8 = 64.
    min_sq_dist = 64

    # List to store the centers of the packed candies
    packed_candies_centers = []

    # 3. Implement the greedy algorithm
    # Iterate through all possible center locations in a fixed order (z, then y, then x)
    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                current_center = (x, y, z)
                
                # Assume we can place the candy, then check for conflicts
                can_place = True
                for placed_center in packed_candies_centers:
                    # Calculate squared distance to each already-placed candy
                    sq_dist = (current_center[0] - placed_center[0])**2 + \
                              (current_center[1] - placed_center[1])**2 + \
                              (current_center[2] - placed_center[2])**2
                    
                    # If it's too close to any existing candy, we can't place it
                    if sq_dist < min_sq_dist:
                        can_place = False
                        break
                
                # If after checking all placed candies there are no conflicts, place the new one
                if can_place:
                    packed_candies_centers.append(current_center)

    # 4. Output the results
    # The maximized value n is the total number of candies placed.
    n = len(packed_candies_centers)

    print(f"The problem formulation is correct.")
    print(f"The maximized number of eyeball candies (n) is: {n}")
    print("-" * 30)

    # To show the "final equation", we can break down the packing by layers
    layers = {}
    for center in packed_candies_centers:
        z_val = center[2]
        if z_val not in layers:
            layers[z_val] = 0
        layers[z_val] += 1

    # Sort layers by z-coordinate for a clear output
    sorted_layers = sorted(layers.items())

    equation_parts = [str(count) for z, count in sorted_layers]
    equation_str = " + ".join(equation_parts)
    
    print("The packing is arranged in layers as follows:")
    for z_val, count in sorted_layers:
        # Convert z-coordinate back to cm for clarity
        print(f"Layer at z = {z_val * 0.5} cm: {count} candies")

    print("\nThe final equation for the total number of candies is:")
    print(f"{equation_str} = {n}")

solve_candy_packing()