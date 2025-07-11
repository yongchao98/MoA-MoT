import random
import math

def solve():
    """
    Calculates the probability that a circle of radius R=6 intersects exactly N=47
    cells of a unit grid, using a Monte Carlo simulation.
    """
    
    R = 6.0
    R_sq = R**2
    target_intersections = 47
    num_samples = 2 * 10**6  # Number of random throws

    print(f"Parameters of the problem:")
    print(f"Circle radius R = {R}")
    print(f"Squared radius R^2 = {R_sq}")
    print(f"Target number of intersected cells N = {target_intersections}")
    print(f"Running Monte Carlo simulation with {num_samples} samples...\n")

    successful_throws = 0

    # The center of the circle (x0, y0) is randomly chosen from a unit square [0, 1) x [0, 1).
    for _ in range(num_samples):
        x0 = random.random()
        y0 = random.random()

        intersected_cells_count = 0
        
        # We only need to check cells that could possibly be intersected.
        # A circle with center in [0,1)x[0,1) and R=6 spans x-coordinates in (-6, 7)
        # and y-coordinates in (-6, 7).
        # A cell (k,m) is the square [k, k+1] x [m, m+1].
        # The indices k and m can range from -6 to 6.
        for k in range(-6, 7):
            for m in range(-6, 7):
                
                # A cell is intersected if it's neither fully inside nor fully outside the circle.
                # This is equivalent to: min_dist <= R <= max_dist from the circle's center
                # to the closest/farthest point of the cell.
                
                # Find the coordinates of the point in the cell [k, k+1]x[m, m+1]
                # that is closest to the circle's center (x0, y0).
                if k >= 1:
                    px_c = k
                elif k == 0:
                    px_c = x0
                else:  # k <= -1
                    px_c = k + 1

                if m >= 1:
                    py_c = m
                elif m == 0:
                    py_c = y0
                else:  # m <= -1
                    py_c = m + 1
                
                # Squared distance to the closest point in the cell.
                dist_sq_min = (px_c - x0)**2 + (py_c - y0)**2

                # If the closest point is farther than R, the cell is entirely outside.
                if dist_sq_min > R_sq:
                    continue

                # Find the farthest point in the cell from the center (x0, y0).
                # This will be one of the four corners of the cell.
                # The squared distance to the 4 corners:
                dist_sq_c1 = (k - x0)**2 + (m - y0)**2
                dist_sq_c2 = (k + 1 - x0)**2 + (m - y0)**2
                dist_sq_c3 = (k - x0)**2 + (m + 1 - y0)**2
                dist_sq_c4 = (k + 1 - x0)**2 + (m + 1 - y0)**2
                
                dist_sq_max = max(dist_sq_c1, dist_sq_c2, dist_sq_c3, dist_sq_c4)
                
                # If the farthest point is closer than R, the cell is entirely inside.
                if dist_sq_max < R_sq:
                    continue
                
                # Otherwise, the circumference intersects the cell.
                intersected_cells_count += 1
        
        if intersected_cells_count == target_intersections:
            successful_throws += 1
            
    probability = successful_throws / num_samples

    print("Simulation Results:")
    print(f"Number of throws where N=47: {successful_throws}")
    print(f"Total number of throws: {num_samples}")
    print(f"The final probability is the ratio of these two numbers.")
    print(f"Equation: P(N=47) = {successful_throws} / {num_samples}")
    print(f"Approximate probability: {probability:.4g}")

solve()