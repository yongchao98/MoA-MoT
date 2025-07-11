import random
import math

def solve_circumference_probability():
    """
    Calculates the probability that a randomly thrown circle of radius 6
    on a unit grid intersects exactly 47 cells using a Monte Carlo simulation.
    """
    R = 6.0
    R2 = R * R
    num_samples = 2000000  # Number of random throws. Larger is more accurate but slower.
    target_intersections = 47
    count_for_target = 0

    # Run the simulation for each random sample
    for _ in range(num_samples):
        # Generate a random center for the circle in the unit square [0,1) x [0,1)
        x = random.random()
        y = random.random()

        n_intersected = 0

        # Define the range of grid cells to check.
        # A cell [k,k+1]x[l,l+1] can only be intersected if it's near the circle.
        k_min = math.floor(x - R)
        k_max = math.ceil(x + R)
        l_min = math.floor(y - R)
        l_max = math.ceil(y + R)

        for k in range(k_min, k_max):
            for l in range(l_min, l_max):
                # The current cell is the square [k, k+1] x [l, l+1]

                # --- Check for intersection ---

                # 1. Calculate the squared distance from the center (x,y) to the
                #    closest point in the current cell.
                # If x is within the cell's x-span [k, k+1], the closest x-distance is 0.
                # Otherwise, it's the distance to the nearest edge (k or k+1).
                dx_min = 0
                if x < k:
                    dx_min = k - x
                elif x > k + 1:
                    dx_min = x - (k + 1)

                dy_min = 0
                if y < l:
                    dy_min = l - y
                elif y > l + 1:
                    dy_min = y - (l + 1)

                min_dist_sq = dx_min**2 + dy_min**2

                # Optimization: If the minimum distance is greater than R, the circle
                # cannot touch the cell, so we can skip to the next cell.
                if min_dist_sq > R2:
                    continue

                # 2. Calculate the squared distance from the center (x,y) to the
                #    farthest point in the cell. This will always be one of the four corners.
                dx1 = k - x
                dx2 = (k + 1) - x
                dy1 = l - y
                dy2 = (l + 1) - y

                max_dist_sq = max(dx1**2 + dy1**2, dx1**2 + dy2**2, dx2**2 + dy1**2, dx2**2 + dy2**2)
                
                # Optimization: If the maximum distance is less than R, the cell is
                # entirely inside the circle, so the circumference doesn't cross it.
                if max_dist_sq < R2:
                    continue
                
                # 3. If R is between min_dist and max_dist, the circumference intersects the cell.
                # (min_dist_sq <= R2 <= max_dist_sq)
                n_intersected += 1

        # Check if the number of intersections for this throw matches our target
        if n_intersected == target_intersections:
            count_for_target += 1

    # Calculate the final probability
    probability = count_for_target / num_samples
    
    # Print the numbers used in the final equation
    print(f"Number of events with {target_intersections} intersections: {count_for_target}")
    print(f"Total number of random throws: {num_samples}")
    print(f"Final equation: {count_for_target} / {num_samples} = {probability}")
    
    # Print the approximate solution with four significant digits
    print(f"\nApproximate probability: {probability:.4f}")
    
if __name__ == '__main__':
    solve_circumference_probability()
