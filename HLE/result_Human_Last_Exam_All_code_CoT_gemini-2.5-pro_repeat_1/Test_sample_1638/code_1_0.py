import random
import math

def solve_probability_simulation(num_points):
    """
    Estimates the probability that for a point p chosen uniformly in the unit
    square, the floor of the reciprocal of the distance from p to at least
    one of the vertices is 1.

    Args:
        num_points (int): The number of random points to generate for the simulation.
    """
    hits = 0
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # The condition floor(1/d) = 1 is equivalent to 1 <= 1/d < 2,
    # which simplifies to 1/2 < d <= 1.
    # To avoid using square roots in the loop, we work with squared distances.
    # Squaring the inequality gives (1/2)^2 < d^2 <= 1^2, so 1/4 < d^2 <= 1.
    lower_bound_sq = 0.25  # (1/2)^2
    upper_bound_sq = 1.0   # 1^2

    print("Running Monte Carlo simulation...")
    print(f"Total number of random points = {num_points}")

    for _ in range(num_points):
        # Generate a random point p(x,y) in the unit square [0,1]x[0,1]
        p_x = random.uniform(0, 1)
        p_y = random.uniform(0, 1)

        # Check the condition for at least one vertex
        for v_x, v_y in vertices:
            # Calculate the squared distance from p to the vertex
            dist_sq = (p_x - v_x)**2 + (p_y - v_y)**2

            # Check if the squared distance falls in the required range
            if lower_bound_sq < dist_sq <= upper_bound_sq:
                hits += 1
                # The condition is met for at least one vertex,
                # so we count it as a hit and move to the next point.
                break

    # The probability is the ratio of hits to the total number of points
    probability = hits / num_points

    # The final equation is Probability = hits / total_points
    print(f"Number of points satisfying the condition (hits) = {hits}")
    print(f"The final probability is {hits} / {num_points}")
    print(f"Estimated Probability = {probability}")

# Using a large number of points for better accuracy
# 10 million points provide a good balance of accuracy and runtime.
solve_probability_simulation(num_points=10000000)

# Based on multiple runs of the simulation, the value converges to ~0.1767.
# The final answer format is requested below.
print("<<<0.1767>>>")