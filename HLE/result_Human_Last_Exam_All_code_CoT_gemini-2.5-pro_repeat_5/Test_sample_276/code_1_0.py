import math
import random

def count_intersected_cells(x, y, R=6):
    """
    Counts the number of grid cells intersected by a circumference.
    
    A cell [n, n+1] x [m, m+1] is intersected if the circumference passes
    through it. This is true if the radius R is between the minimum and
    maximum distance from the circle's center (x,y) to any point in the cell.
    """
    R_sq = R * R
    count = 0
    
    # Determine the range of integer cell coordinates to check
    n_min = math.floor(x - R)
    n_max = math.ceil(x + R)
    m_min = math.floor(y - R)
    m_max = math.ceil(y + R)

    for n in range(n_min, n_max):
        for m in range(m_min, m_max):
            # Cell is the square [n, n+1] x [m, m+1]

            # Find the squared distance from (x,y) to the closest point in the cell
            closest_x = max(n, min(x, n + 1))
            closest_y = max(m, min(y, m + 1))
            min_dist_sq = (x - closest_x)**2 + (y - closest_y)**2

            # Find the squared distance from (x,y) to the farthest corner of the cell
            # This is one of the four corners of the cell square.
            d1 = (x - n)**2 + (y - m)**2
            d2 = (x - (n + 1))**2 + (y - m)**2
            d3 = (x - n)**2 + (y - (m + 1))**2
            d4 = (x - (n + 1))**2 + (y - (m + 1))**2
            max_dist_sq = max(d1, d2, d3, d4)

            if min_dist_sq <= R_sq and R_sq <= max_dist_sq:
                count += 1
    return count

def run_simulation(num_trials=100000):
    """
    Runs the Monte Carlo simulation to find the probability.
    """
    R = 6
    target_intersections = 47
    count_for_target = 0

    for _ in range(num_trials):
        # Generate a random center (x,y) in the unit square [0,1]x[0,1]
        x_center = random.uniform(0, 1)
        y_center = random.uniform(0, 1)

        num_intersected = count_intersected_cells(x_center, y_center, R)
        
        if num_intersected == target_intersections:
            count_for_target += 1
            
    probability = count_for_target / num_trials
    
    print(f"Number of trials: {num_trials}")
    print(f"Number of times 47 cells were intersected: {count_for_target}")
    print(f"The equation for the probability is: P(K=47) = {count_for_target} / {num_trials}")
    print(f"Approximate probability: {probability:.4f}")

if __name__ == '__main__':
    # A smaller number of trials for quick execution, but a larger number
    # (e.g., 100,000 or more) is needed for a more reliable estimate.
    # Based on mathematical properties, the number of intersected cells is
    # always even for a circle of radius 6. Therefore, the count for 47
    # intersections is expected to be 0.
    run_simulation(num_trials=10000)
