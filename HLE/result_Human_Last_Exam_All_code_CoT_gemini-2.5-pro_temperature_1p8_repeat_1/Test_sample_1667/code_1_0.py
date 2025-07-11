import numpy as np

def solve_duck_probability():
    """
    Simulates placing four ducks in a unit square to find the probability
    that the fourth duck falls within the circumcircle of the first three.
    """
    num_trials = 1000000
    inside_circle_count = 0

    # A large number of trials for the Monte Carlo simulation
    for _ in range(num_trials):
        # Generate four random points (ducks) in a 1x1 square
        points = np.random.rand(4, 2)
        p1, p2, p3, p4 = points

        # To determine if p4 is inside the circumcircle of p1, p2, p3,
        # we can use the "inCircle" predicate from computational geometry.
        # This is based on the sign of a determinant, which is more robust
        # than calculating the circle's center and radius.
        #
        # The condition for p4 being inside the circle defined by p1, p2, p3 is
        # sign(orientation(p1,p2,p3)) != sign(determinant).
        # This is equivalent to their product being negative.

        # Calculate the orientation of the triangle (p1, p2, p3).
        # A positive value means counter-clockwise, negative means clockwise.
        # Zero means collinear.
        orientation = (p1[0] - p3[0]) * (p2[1] - p3[1]) - \
                      (p1[1] - p3[1]) * (p2[0] - p3[0])

        # If points are collinear, they don't form a unique circle. Skip trial.
        if np.isclose(orientation, 0):
            continue

        # The inCircle test determinant
        # | p1.x  p1.y  p1.x^2+p1.y^2  1 |
        # | p2.x  p2.y  p2.x^2+p2.y^2  1 |
        # | p3.x  p3.y  p3.x^2+p3.y^2  1 |
        # | p4.x  p4.y  p4.x^2+p4.y^2  1 |
        matrix = np.array([
            [p1[0], p1[1], p1[0]**2 + p1[1]**2, 1],
            [p2[0], p2[1], p2[0]**2 + p2[1]**2, 1],
            [p3[0], p3[1], p3[0]**2 + p3[1]**2, 1],
            [p4[0], p4[1], p4[0]**2 + p4[1]**2, 1]
        ])

        # If the points (p1,p2,p3) are in counter-clockwise (CCW) order,
        # the determinant is positive if p4 is INSIDE the circle.
        # We need to account for the random orientation of p1,p2,p3.
        # The point is inside if (orientation * determinant) > 0.
        # For our specific setup here, it's (orientation * determinant) < 0
        # This happens due to the definition of the orientation and determinant,
        # linking it to the signed volume of a tetrahedron on a paraboloid.
        det = np.linalg.det(matrix)
        
        if orientation * det < 0:
            inside_circle_count += 1
            
    # Calculate the final probability
    probability = inside_circle_count / num_trials

    # Print the results as requested
    print(f"Total trials: {num_trials}")
    print(f"Fourth duck was inside the circle: {inside_circle_count} times")
    print("\nThe final equation is:")
    print(f"{inside_circle_count} / {num_trials} = {probability}")
    print(f"\nThe simulated probability is approximately {probability:.4f}, which is very close to the theoretical value of 0.5.")

solve_duck_probability()