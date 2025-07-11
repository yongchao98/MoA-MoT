import random
import math

def count_intersected_cells(x, y, r=6):
    """
    Counts the number of grid cells intersected by a circumference.
    
    Args:
        x (float): The x-coordinate of the circle's center, in [0, 1].
        y (float): The y-coordinate of the circle's center, in [0, 1].
        r (float): The radius of the circle.
        
    Returns:
        int: The number of intersected grid cells.
    """
    count = 0
    r_sq = r**2
    
    # Define a bounding box of grid cells to check. Since the center (x,y)
    # is in [0,1]x[0,1] and r=6, the circle can span from -6 to 7.
    # We add a small margin for safety.
    k_min = math.floor(x - r) - 1
    k_max = math.ceil(x + r) + 1
    l_min = math.floor(y - r) - 1
    l_max = math.ceil(y + r) + 1

    for k in range(k_min, k_max):
        for l in range(l_min, l_max):
            # The current grid cell is the square [k, k+1] x [l, l+1].
            
            # Find the closest point in the cell to the circle's center (x,y).
            # This is done by clamping the coordinates of the center to the cell's boundaries.
            px_min = max(k, min(x, k + 1))
            py_min = max(l, min(y, l + 1))
            d_min_sq = (px_min - x)**2 + (py_min - y)**2
            
            # If the squared minimum distance is greater than r^2, the circle
            # is entirely outside the cell, so it doesn't intersect.
            if d_min_sq > r_sq:
                continue
            
            # Find the furthest point in the cell to (x,y). This will be one of the four corners.
            corners = [(k, l), (k + 1, l), (k, l + 1), (k + 1, l + 1)]
            d_max_sq = 0
            for cx, cy in corners:
                d_sq = (cx - x)**2 + (cy - y)**2
                if d_sq > d_max_sq:
                    d_max_sq = d_sq

            # If the squared maximum distance is less than r^2, the cell is
            # entirely inside the circle, so the circumference doesn't intersect it.
            if d_max_sq < r_sq:
                continue

            # If neither of the above, the circumference must intersect the cell.
            count += 1
            
    return count

def main():
    """
    Runs the Monte Carlo simulation and prints the result.
    """
    # A large number of samples for a more accurate result.
    num_samples = 2000000 
    
    # Count of trials where the number of intersections is 47.
    favorable_outcomes = 0
    
    # Use a fixed seed for reproducibility of the result.
    random.seed(0)
    
    for _ in range(num_samples):
        # Generate a random center for the circle in the [0, 1] x [0, 1] square.
        rand_x = random.random()
        rand_y = random.random()
        
        # Count the number of intersected cells.
        num_intersections = count_intersected_cells(rand_x, rand_y)
        
        # Check if the outcome is the one we are interested in.
        if num_intersections == 47:
            favorable_outcomes += 1
            
    # Calculate the probability.
    probability = favorable_outcomes / num_samples
    
    print("This script calculates the probability using a Monte Carlo simulation.")
    print(f"Number of favorable outcomes (intersections = 47): {favorable_outcomes}")
    print(f"Total number of samples: {num_samples}")
    print(f"The final equation is: {favorable_outcomes} / {num_samples} = {probability}")
    print(f"The approximate probability with four significant digits is {probability:.4g}.")

if __name__ == "__main__":
    main()
