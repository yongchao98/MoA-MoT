import math
import random

def count_intersected_cells(x, y, R):
    """
    Counts the number of grid cells intersected by a circle of radius R
    centered at (x, y), where (x,y) is within the [0,1]x[0,1] square.
    """
    R_sq = R * R
    count = 0
    
    # Determine the range of grid cells to check. A cell S_km is outside
    # the circle's reach if its minimum distance to (x,y) is > R.
    # The minimum distance from (x,y) in [0,1]^2 to S_km=[k,k+1]x[m,m+1] grows with |k| and |m|.
    # We add a small buffer to the range for safety.
    max_coord = int(R) + 2
    
    for k in range(-max_coord, max_coord):
        for m in range(-max_coord, max_coord):
            # The grid cell is S_km = [k, k+1] x [m, m+1].
            
            # Find the point (px, py) in the cell S_km closest to (x, y).
            px = max(k, min(k + 1, x))
            py = max(m, min(m + 1, y))
            
            min_dist_sq = (px - x)**2 + (py - y)**2
            
            # Check if the minimum distance is less than or equal to R.
            # If not, the circle is too far to intersect this cell.
            if min_dist_sq <= R_sq:
                # Find the point in the cell farthest from (x, y). This must be one of the four corners.
                farthest_x = k if abs(k - x) > abs(k + 1 - x) else k + 1
                farthest_y = m if abs(m - y) > abs(m + 1 - y) else m + 1
                
                max_dist_sq = (farthest_x - x)**2 + (farthest_y - y)**2

                # To intersect, the circle must not completely contain or be contained by the cell.
                # Since R=6 and cell is 1x1, R is large enough that the circumference
                # must pass through if the cell is close enough. The condition is
                # that R is between the min and max distance to the cell.
                if max_dist_sq >= R_sq:
                    count += 1
    return count

def main():
    """
    Runs a Monte Carlo simulation to find the probability.
    """
    R = 6.0
    TARGET_CELL_COUNT = 47
    
    # Using a large number of samples for better precision.
    num_samples = 5_000_000
    
    # Count how many times the number of intersected cells is 47
    favorable_outcomes = 0
    
    # Due to grid symmetry, sampling from a quadrant [0, 0.5]x[0, 0.5] is sufficient
    # and equivalent to sampling the whole [0, 1]x[0, 1] square.
    for _ in range(num_samples):
        # Generate a random center point for the circle.
        rand_x = random.uniform(0, 0.5)
        rand_y = random.uniform(0, 0.5)
        
        num_cells = count_intersected_cells(rand_x, rand_y, R)
        
        if num_cells == TARGET_CELL_COUNT:
            favorable_outcomes += 1
            
    # Calculate the probability
    probability = favorable_outcomes / num_samples
    
    # The result of the simulation might be exactly zero.
    # In that case, we can infer the true probability is likely zero.
    print(f"Based on {num_samples} samples, the estimated probability is:")
    # Print with enough precision to see small non-zero results.
    print(f"P(N={TARGET_CELL_COUNT}) = {probability:.8f}")

    # Further analysis of boundary conditions
    # We can check the number of cells for a center on the boundary explicitly.
    n_tangent = count_intersected_cells(0, 0.5, R)
    print(f"\nFor a center on a tangency line, e.g., (0, 0.5), the number of intersected cells is: {n_tangent}")
    
    # The approximate solution as requested. The simulation will likely yield 0.
    # If a non-zero answer is expected, it would be extremely small. Given the problem's nature, 0 is the most likely rigorous answer.
    # We will format the simulation result to 4 significant digits.
    approx_solution = f"{probability:.4g}"
    print(f"\nApproximate solution with four significant digits: {approx_solution}")


if __name__ == '__main__':
    main()
