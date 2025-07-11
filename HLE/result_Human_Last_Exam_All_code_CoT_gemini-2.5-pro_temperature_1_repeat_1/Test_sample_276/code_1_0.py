import math
import random

def count_intersections(cx, cy, R):
    """
    Counts the number of grid cells intersected by the circumference of a circle
    with a given center (cx, cy) and radius R.
    """
    count = 0
    # Define a bounding box of cells to check, with a small margin
    m_min = math.floor(cx - R) - 1
    m_max = math.ceil(cx + R) + 1
    n_min = math.floor(cy - R) - 1
    n_max = math.ceil(cy + R) + 1

    R_squared = R**2

    for m in range(m_min, m_max):
        for n in range(n_min, n_max):
            # Check for intersection with the grid cell [m, m+1] x [n, n+1]
            
            # 1. Find the squared distance from the circle's center to the closest point in the cell.
            # This determines if the circle's disk intersects the cell at all.
            dx_min = 0
            if cx < m:
                dx_min = m - cx
            elif cx > m + 1:
                dx_min = cx - (m + 1)
            
            dy_min = 0
            if cy < n:
                dy_min = n - cy
            elif cy > n + 1:
                dy_min = cy - (n + 1)

            dist_sq_min = dx_min**2 + dy_min**2
            
            # If the closest point of the cell is farther than R, the cell is entirely outside the circle.
            if dist_sq_min > R_squared:
                continue

            # 2. Find the squared distance to the farthest point in the cell.
            # This determines if the cell is completely inside the circle.
            dx_max = max(abs(cx - m), abs(cx - (m + 1)))
            dy_max = max(abs(cy - n), abs(cy - (n + 1)))
            dist_sq_max = dx_max**2 + dy_max**2
            
            # If the farthest point is closer than R, the cell is entirely inside the circle's disk,
            # so the circumference does not intersect it.
            if dist_sq_max < R_squared:
                continue
            
            # If the cell is neither completely outside nor completely inside,
            # the circumference must intersect it.
            count += 1
            
    return count

def main():
    """
    Main function to run the Monte Carlo simulation.
    """
    R = 6
    target_intersections = 47
    # A large number of trials is needed for precision to four significant digits.
    num_trials = 2000000

    hits = 0
    # Run the simulation
    for _ in range(num_trials):
        # Generate a random center (cx, cy) in the unit square [0, 1] x [0, 1]
        cx = random.random()
        cy = random.random()
        
        intersections = count_intersections(cx, cy, R)
        if intersections == target_intersections:
            hits += 1

    # Calculate and print the final probability
    probability = hits / num_trials
    
    print(f"Radius (R): {R}")
    print(f"Target number of intersections: {target_intersections}")
    print(f"Total trials: {num_trials}")
    print(f"Number of times the circumference intersected {target_intersections} cells: {hits}")
    print("\nThe probability is the ratio of hits to total trials:")
    print(f"P(N={target_intersections}) = {hits} / {num_trials} = {probability:.4f}")

if __name__ == "__main__":
    main()