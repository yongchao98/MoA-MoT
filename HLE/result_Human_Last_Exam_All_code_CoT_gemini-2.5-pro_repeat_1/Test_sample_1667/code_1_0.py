import numpy as np

def solve_duck_problem():
    """
    This function calculates the probability that a fourth randomly placed duck
    falls within the circumcircle of three other randomly placed ducks in a unit square.
    It uses a Monte Carlo simulation to estimate the probability.
    """
    num_trials = 1000000
    inside_circle_count = 0

    # Generate all random points at once for efficiency
    # 4 points per trial, 2 coordinates (x,y) each
    points = np.random.rand(num_trials, 4, 2)

    for i in range(num_trials):
        p1, p2, p3, p4 = points[i]

        # To robustly check if p4 is in the circumcircle of p1, p2, p3,
        # we can use a determinant test. The sign of this determinant
        # relates to the orientation of the 4 points lifted to a 3D paraboloid.
        
        # Matrix for the orientation test of p1, p2, p3
        # This determines if the triangle is wound clockwise or counter-clockwise.
        # sign(det(M_orient)) tells us the orientation.
        M_orient = np.array([
            [p1[0], p1[1], 1],
            [p2[0], p2[1], 1],
            [p3[0], p3[1], 1]
        ])
        
        # Matrix for the in-circle test.
        # sign(det(M_incircle)) tells if p4 is inside or outside the circle,
        # but its interpretation depends on the orientation of p1,p2,p3.
        M_incircle = np.array([
            [p1[0], p1[1], p1[0]**2 + p1[1]**2, 1],
            [p2[0], p2[1], p2[0]**2 + p2[1]**2, 1],
            [p3[0], p3[1], p3[0]**2 + p3[1]**2, 1],
            [p4[0], p4[1], p4[0]**2 + p4[1]**2, 1]
        ])

        try:
            det_orient = np.linalg.det(M_orient)
            det_incircle = np.linalg.det(M_incircle)
        except np.linalg.LinAlgError:
            # This is extremely unlikely with random float coordinates
            continue

        # If p1,p2,p3 are collinear, det_orient is 0. We skip these trials.
        if det_orient == 0:
            continue

        # The condition for p4 to be inside the circumcircle is that the sign of
        # the orientation determinant and the incircle determinant match.
        # (Technically, `det_incircle` has a sign convention opposite to `det_orient`,
        # so we check if their product is negative).
        # Let's check a known configuration: p1=(0,0), p2=(2,0), p3=(1,1), p4=(1,0)
        # orient: CCW (+)
        # incircle: center (1,-..), radius sqrt(2). (1,0) is inside.
        # det_incircle must be positive for CCW. so product > 0.
        if det_orient * det_incircle > 0:
            inside_circle_count += 1
            
    probability = inside_circle_count / num_trials

    print(f"Based on a simulation with {num_trials} trials:")
    print(f"Estimated Probability: {probability:.6f}")
    print("\nThe final probability is estimated as the ratio of successful trials to the total number of trials.")
    print(f"Number of successful trials (duck #4 in the circle): {inside_circle_count}")
    print(f"Total number of trials: {num_trials}")
    print(f"Final equation: {inside_circle_count} / {num_trials} = {probability:.6f}")

solve_duck_problem()