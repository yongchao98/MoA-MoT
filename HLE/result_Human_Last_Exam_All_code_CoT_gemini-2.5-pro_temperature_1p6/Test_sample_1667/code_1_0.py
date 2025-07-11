import numpy as np

def solve_probability():
    """
    Runs a Monte Carlo simulation to find the probability that a fourth random
    point in a unit square falls within the circumcircle of the first three.
    """
    num_trials = 1_000_000
    inside_circle_count = 0

    # Generate all random points at once for efficiency
    points = np.random.rand(num_trials, 4, 2)

    for i in range(num_trials):
        p1, p2, p3, p4 = points[i]

        # Use the linear algebra method to find the circumcenter (h, k)
        # of the triangle formed by p1, p2, p3.
        # The system of equations for the center (h,k) is:
        # 2(x2-x1)h + 2(y2-y1)k = x2^2+y2^2 - x1^2-y1^2
        # 2(x3-x2)h + 2(y3-y2)k = x3^2+y3^2 - x2^2-y2^2
        
        # To avoid verbose indexing, unpack coordinates
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3
        x4, y4 = p4
        
        # Coefficients of the linear system
        a = 2 * (x2 - x1)
        b = 2 * (y2 - y1)
        c = x2**2 + y2**2 - x1**2 - y1**2
        
        d = 2 * (x3 - x2)
        e = 2 * (y3 - y2)
        f = x3**2 + y3**2 - x2**2 - y2**2
        
        # Solve for h, k using Cramer's rule. The denominator is the determinant.
        determinant = a * e - b * d
        
        # If the determinant is near zero, the points are collinear.
        # This is extremely rare with random floats, so we skip the trial.
        if abs(determinant) < 1e-10:
            continue
            
        h = (c * e - b * f) / determinant
        k = (a * f - c * d) / determinant
        
        # Calculate radius squared from the center to any of the three points (e.g., p1)
        radius_sq = (x1 - h)**2 + (y1 - k)**2
        
        # Calculate the distance squared from the fourth point to the center
        dist_sq = (x4 - h)**2 + (y4 - k)**2
        
        # Check if the fourth point is inside the circle
        if dist_sq < radius_sq:
            inside_circle_count += 1
            
    # Calculate the final probability
    probability = inside_circle_count / num_trials
    
    # Per the instructions, output the final equation
    print(f"The simulation involved {num_trials} trials.")
    print(f"The fourth duck was inside the circle {inside_circle_count} times.")
    print(f"The calculated probability is:")
    print(f"{inside_circle_count} / {num_trials} = {probability}")

solve_probability()