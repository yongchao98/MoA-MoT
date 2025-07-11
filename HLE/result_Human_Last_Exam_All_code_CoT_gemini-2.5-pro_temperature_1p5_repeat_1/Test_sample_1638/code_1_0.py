import random
import math

def solve_probability_simulation():
    """
    This function calculates the probability that for a point p, chosen
    uniformly at random in the unit square, the floor of the reciprocal of
    the distance from p to at least one of the vertices is 1.

    This is done using a Monte Carlo simulation.
    """
    num_trials = 10000000  # Number of random points to generate for accuracy
    success_count = 0      # Counter for points satisfying the condition

    # The four vertices of the unit square
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # Loop for each trial
    for _ in range(num_trials):
        # Generate a random point p(x, y) in the unit square [0,1] x [0,1]
        p_x = random.random()
        p_y = random.random()

        # Check the condition for each of the four vertices
        for v_x, v_y in vertices:
            # Calculate the squared distance from the point to the vertex
            # The condition is 1/2 < d <= 1, which is equivalent to 1/4 < d^2 <= 1
            # Using squared distance is more efficient as it avoids the sqrt operation.
            dist_sq = (p_x - v_x)**2 + (p_y - v_y)**2

            if 1/4 < dist_sq <= 1:
                # The condition is met for at least one vertex.
                # Increment count and break the inner loop to check the next point.
                success_count += 1
                break

    # The probability is the ratio of successful trials to the total number of trials
    probability = success_count / num_trials

    # Output the components of the final calculation
    print(f"Total points sampled (num_trials): {num_trials}")
    print(f"Points satisfying the condition (success_count): {success_count}")
    print(f"The final probability is the result of the equation:")
    print(f"P = success_count / num_trials")
    print(f"P = {success_count} / {num_trials}")
    print(f"P ≈ {probability}")


solve_probability_simulation()
<<<0.766>>>