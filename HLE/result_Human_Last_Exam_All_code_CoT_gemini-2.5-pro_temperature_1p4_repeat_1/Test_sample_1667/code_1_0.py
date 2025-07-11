import numpy as np

def solve():
    """
    Runs a Monte Carlo simulation to find the probability that a fourth random
    point in a unit square falls inside the circumcircle of three other
    random points.
    """
    num_trials = 100000
    inside_circle_count = 0

    for _ in range(num_trials):
        # Generate four random points in the unit square
        points = np.random.rand(4, 2)
        p1, p2, p3, p4 = points

        # Unpack coordinates for clarity
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3
        x4, y4 = p4

        # Calculate the denominator of the circumcenter formula.
        # This is non-zero if the points are not collinear.
        D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

        if abs(D) < 1e-10:  # Points are nearly collinear, skip this trial
            continue

        # Calculate the coordinates of the circumcenter (cx, cy)
        s1, s2, s3 = x1**2 + y1**2, x2**2 + y2**2, x3**2 + y3**2
        cx = (s1 * (y2 - y3) + s2 * (y3 - y1) + s3 * (y1 - y2)) / D
        cy = (s1 * (x3 - x2) + s2 * (x1 - x3) + s3 * (x2 - x1)) / D

        # Calculate the squared radius of the circumcircle
        radius_sq = (x1 - cx)**2 + (y1 - cy)**2

        # Calculate the squared distance of the fourth point from the circumcenter
        dist_sq = (x4 - cx)**2 + (y4 - cy)**2

        # Check if the fourth point is inside the circle
        if dist_sq < radius_sq:
            inside_circle_count += 1
            
    # Calculate the estimated probability
    probability = inside_circle_count / num_trials

    print(f"Number of trials: {num_trials}")
    print(f"Number of times the 4th duck was inside the circle: {inside_circle_count}")
    # This line prints the numbers in the final calculation as requested
    print(f"Final calculation: {inside_circle_count} / {num_trials} = {probability}")
    print("The theoretical probability is 0.5")

solve()
