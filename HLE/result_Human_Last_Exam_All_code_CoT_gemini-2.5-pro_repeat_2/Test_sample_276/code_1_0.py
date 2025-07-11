import math
import random

def count_intersections(x_c, y_c, R):
    """
    Counts the number of grid cells intersected by a circle with center (x_c, y_c) and radius R.
    """
    R_sq = R * R
    count = 0
    
    # Determine the range of grid cells to check based on the circle's bounding box.
    n_min = math.floor(x_c - R)
    n_max = math.ceil(x_c + R)
    m_min = math.floor(y_c - R)
    m_max = math.ceil(y_c + R)

    # Iterate over each potentially intersected cell [n, n+1] x [m, m+1].
    for n in range(n_min, n_max):
        for m in range(m_min, m_max):
            
            # Find the squared distance from the center (x_c, y_c) to the closest point in the cell.
            closest_x = max(n, min(x_c, n + 1))
            closest_y = max(m, min(y_c, m + 1))
            d_min_sq = (x_c - closest_x)**2 + (y_c - closest_y)**2

            # If the circle is entirely outside the cell, skip to the next cell.
            if d_min_sq >= R_sq:
                continue

            # Find the squared distance from the center to the furthest point in the cell (one of the four corners).
            d_max_sq = max(
                (n - x_c)**2 + (m - y_c)**2,
                ((n + 1) - x_c)**2 + (m - y_c)**2,
                (n - x_c)**2 + ((m + 1) - y_c)**2,
                ((n + 1) - x_c)**2 + ((m + 1) - y_c)**2
            )

            # The circumference intersects the cell's interior if R is strictly between the min and max distance.
            if R_sq < d_max_sq:
                count += 1
                
    return count

def solve_geometric_probability():
    """
    Runs a Monte Carlo simulation to find the probability.
    """
    R = 6
    TARGET_N = 47
    
    # A large number of trials is needed for precision as the probability is small.
    # This may take a minute or two to execute.
    NUM_TRIALS = 5000000 

    success_count = 0
    
    # Use a fixed seed for reproducibility.
    random.seed(42)
    
    for i in range(NUM_TRIALS):
        # Generate a random center (x_c, y_c) in the unit square [0, 1] x [0, 1].
        x_c = random.random()
        y_c = random.random()

        # Count the number of intersected cells.
        num_intersections = count_intersections(x_c, y_c, R)
        
        # Check if the count matches the target.
        if num_intersections == TARGET_N:
            success_count += 1
    
    # The probability is the ratio of successful trials to the total number of trials.
    probability = success_count / NUM_TRIALS
    
    # Output the final equation and the result.
    print(f"Radius (R): {R}")
    print(f"Target number of intersections (N): {TARGET_N}")
    print("\n--- Simulation Results ---")
    print(f"Total trials: {NUM_TRIALS}")
    print(f"Successful trials (N = {TARGET_N}): {success_count}")
    print(f"\nProbability = {success_count} / {NUM_TRIALS}")
    print(f"Approximate probability: {probability:.4g}")

# Run the solver.
solve_geometric_probability()
