import random
from collections import Counter

def solve_packing():
    """
    Solves the discrete sphere packing problem using a randomized greedy algorithm.
    """
    # Problem parameters from the user's formulation.
    # The sphere radius is 2cm, and the grid step is 0.5cm.
    # In grid units, the sphere radius is 2 / 0.5 = 4.
    # The minimum distance between centers is the diameter, 2 * 4 = 8 grid units.
    # The squared minimum distance is 8*8 = 64, which matches the user's (4+4)^2.
    min_dist_sq_grid = 64

    # Valid range for sphere centers in grid units, as defined in the problem.
    # x, y in [4, 20], z in [4, 18]
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)

    # Generate all possible center points on the grid
    possible_centers = []
    for x in x_range:
        for y in y_range:
            for z in z_range:
                possible_centers.append((x, y, z))

    # --- Randomized Greedy Search ---
    # We run the search multiple times with random starting orders
    # to find a good solution.
    
    max_n = 0
    best_config = []
    num_trials = 250 # Number of random attempts for a robust result
    
    # Using a fixed seed ensures the result is reproducible.
    random.seed(42)

    for _ in range(num_trials):
        # Shuffle the list of possible centers for a random greedy approach
        random.shuffle(possible_centers)
        
        placed_spheres = []
        for center in possible_centers:
            is_valid = True
            # Check for collision with already placed spheres
            for placed in placed_spheres:
                dist_sq = (center[0] - placed[0])**2 + \
                          (center[1] - placed[1])**2 + \
                          (center[2] - placed[2])**2
                if dist_sq < min_dist_sq_grid:
                    is_valid = False
                    break
            
            if is_valid:
                placed_spheres.append(center)
        
        # If this trial found a better solution, save it
        if len(placed_spheres) > max_n:
            max_n = len(placed_spheres)
            best_config = placed_spheres

    # --- Analyze and Print the Result ---
    # To "output each number in the final equation," we show the composition
    # of the final number by counting spheres on each z-level.

    if not best_config:
        print("Could not place any spheres.")
        return

    # Count spheres at each z-level in the best configuration found
    z_coords = [p[2] for p in best_config]
    z_counts = Counter(z_coords)
    
    # Sort the counts by z-level for a clear, ordered equation
    sorted_counts = sorted(z_counts.items())
    
    # Build the equation string from the counts per layer
    count_terms = [str(count) for z, count in sorted_counts]
    equation_str = " + ".join(count_terms)
    
    # Print the final equation and the total
    print(f"{equation_str} = {max_n}")

# Execute the function to find and print the solution
solve_packing()