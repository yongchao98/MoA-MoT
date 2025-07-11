import math
import random

def count_intersections(x0, y0, r):
    """
    Counts the number of grid cells intersected by a circumference.
    The circumference is centered at (x0, y0) with radius r.
    The grid is a square grid of mesh 1.
    """
    r_squared = r * r
    i_min = math.floor(x0 - r)
    i_max = math.floor(x0 + r)
    j_min = math.floor(y0 - r)
    j_max = math.floor(y0 + r)

    count = 0
    for i in range(i_min, i_max + 1):
        for j in range(j_min, j_max + 1):
            # For cell (i, j), which is the square [i, i+1] x [j, j+1]
            
            # Calculate minimum distance squared from center to cell
            dx = 0
            if x0 < i:
                dx = i - x0
            elif x0 > i + 1:
                dx = x0 - (i + 1)
            
            dy = 0
            if y0 < j:
                dy = j - y0
            elif y0 > j + 1:
                dy = y0 - (j + 1)
            
            dist_sq_min = dx*dx + dy*dy

            # If the center is inside the cell, dist_sq_min is 0.
            # A cell containing the center is always intersected unless r=0.
            if dist_sq_min > r_squared:
                continue

            # Calculate maximum distance squared from center to cell
            corners = [(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)]
            dist_sq_max = 0
            for cx, cy in corners:
                dist_sq = (cx - x0)**2 + (cy - y0)**2
                if dist_sq > dist_sq_max:
                    dist_sq_max = dist_sq
            
            # Check for intersection
            if r_squared <= dist_sq_max:
                count += 1
    return count

def calculate_probability():
    """
    Runs a Monte Carlo simulation to find the probability.
    """
    radius = 6.0
    target_intersections = 47
    
    # A large number of trials is needed for precision.
    # 2,000,000 provides a balance between accuracy and computation time.
    num_trials = 2000000
    
    success_count = 0
    
    # Using a fixed seed for reproducibility of the result.
    random.seed(123)
    
    for _ in range(num_trials):
        x0 = random.random()
        y0 = random.random()
        
        intersections = count_intersections(x0, y0, radius)
        
        if intersections == target_intersections:
            success_count += 1
            
    probability = success_count / num_trials
    
    print(f"The final probability is the result of the equation: {success_count} / {num_trials}")
    print(f"Approximate probability: {probability:.4f}")

calculate_probability()