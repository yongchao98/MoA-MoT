import math
import random

def count_intersections(x0, y0, R):
    """
    Counts the number of grid cells intersected by a circle.

    A cell is intersected by the circumference if it's neither entirely inside
    nor entirely outside the circle.
    
    Args:
        x0 (float): The x-coordinate of the circle's center.
        y0 (float): The y-coordinate of the circle's center.
        R (float): The radius of the circle.

    Returns:
        int: The number of intersected grid cells.
    """
    R_sq = R * R
    intersected_count = 0

    # Determine a bounding box of grid cells to check
    m_min = math.floor(x0 - R)
    m_max = math.ceil(x0 + R)
    n_min = math.floor(y0 - R)
    n_max = math.ceil(y0 + R)

    for m in range(m_min, m_max):
        for n in range(n_min, n_max):
            # The cell is the square [m, m+1] x [n, n+1]
            
            # Find the squared distance from the circle's center to the closest point in the cell.
            # Source: https://stackoverflow.com/questions/5254838/
            dist_x_sq, dist_y_sq = 0, 0
            if x0 < m:
                dist_x_sq = (m - x0)**2
            elif x0 > m + 1:
                dist_x_sq = (x0 - (m + 1))**2
            
            if y0 < n:
                dist_y_sq = (n - y0)**2
            elif y0 > n + 1:
                dist_y_sq = (y0 - (n + 1))**2
                
            min_dist_sq = dist_x_sq + dist_y_sq
            
            # If the closest point in the cell is further than the radius,
            # the entire cell is outside the circle, so it's not intersected.
            if min_dist_sq > R_sq:
                continue
                
            # Find the squared distance to the furthest point in the cell.
            # This will be one of the four corners of the cell.
            max_dist_sq = max((m - x0)**2, ((m + 1) - x0)**2) + max((n - y0)**2, ((n + 1) - y0)**2)

            # If the furthest point in the cell is closer than the radius,
            # the entire cell is inside the circle, so the circumference doesn't intersect it.
            if max_dist_sq < R_sq:
                continue

            # If the cell is not entirely outside and not entirely inside,
            # the circumference must intersect it.
            intersected_count += 1
            
    return intersected_count

def run_simulation(radius, target_intersections, num_samples):
    """
    Runs a Monte Carlo simulation to find the probability of a specific intersection count.
    
    Args:
        radius (float): The radius of the circle.
        target_intersections (int): The target number of intersections.
        num_samples (int): The number of random trials to run.
        
    Returns:
        float: The estimated probability.
    """
    success_count = 0
    
    for _ in range(num_samples):
        # The center of the circle is uniformly distributed.
        # We can analyze its position within a single unit square [0,1]x[0,1].
        x_center = random.random()
        y_center = random.random()
        
        num_intersected = count_intersections(x_center, y_center, radius)
        
        if num_intersected == target_intersections:
            success_count += 1
            
    probability = success_count / num_samples
    return probability

if __name__ == '__main__':
    RADIUS = 6
    TARGET_COUNT = 47
    
    # Using a large number of samples for better accuracy.
    # Note: This can take some time to run.
    NUM_SAMPLES = 2000000 
    
    # Calculate the probability
    prob = run_simulation(RADIUS, TARGET_COUNT, NUM_SAMPLES)
    
    # Print the result with four significant digits
    print(f"The number of cells intersected by the circumference is equal to {TARGET_COUNT}.")
    print(f"Radius of the circumference: {RADIUS}")
    print(f"Number of random samples in simulation: {NUM_SAMPLES}")
    print(f"Estimated probability: {prob:.4g}")
