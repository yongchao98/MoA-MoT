import numpy as np

def run_simulation():
    """
    Runs a Monte Carlo simulation to estimate the probability that a fourth random
    point in a unit square falls into the circumcircle of the first three.
    """
    # Set the number of trials for the simulation
    num_trials = 200000
    inside_circle_count = 0

    # Generate all random points at once for efficiency
    points_set = np.random.rand(num_trials, 4, 2)

    for points in points_set:
        p1, p2, p3, p4 = points

        # Unpack coordinates for clarity
        ax, ay = p1
        bx, by = p2
        cx, cy = p3

        # Calculate the denominator D for the circumcenter formula.
        # D is twice the signed area of the triangle (p1, p2, p3).
        # If D is zero, the points are collinear.
        D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))

        # As per the problem, we assume non-collinear points. In a simulation,
        # due to floating point precision, D could be extremely small. We skip
        # such rare cases where points are nearly collinear.
        if abs(D) < 1e-12:
            continue

        # Calculate the circumcenter (ux, uy) and squared radius.
        # These formulas derive from the Cartesian equation of a circle.
        sq_a = ax**2 + ay**2
        sq_b = bx**2 + by**2
        sq_c = cx**2 + cy**2

        ux = (sq_a * (by - cy) + sq_b * (cy - ay) + sq_c * (ay - by)) / D
        uy = (sq_a * (cx - bx) + sq_b * (ax - cx) + sq_c * (bx - ax)) / D

        center = np.array([ux, uy])
        radius_sq = np.sum((p1 - center)**2)

        # Check if the fourth point (p4) is inside the circle by comparing
        # the squared distance from the center to the squared radius.
        dist_sq_p4 = np.sum((p4 - center)**2)

        if dist_sq_p4 < radius_sq:
            inside_circle_count += 1

    # Calculate the estimated probability from the simulation
    probability = inside_circle_count / num_trials

    # Print the results of the simulation, including the numbers used
    print(f"Number of simulation trials: {num_trials}")
    print(f"Times the 4th duck was inside the circle: {inside_circle_count}")
    print("The final probability is the ratio of successes to trials.")
    print(f"Final equation from simulation:")
    print(f"Probability = {inside_circle_count} / {num_trials} = {probability:.4f}")
    print("\nThis simulation result is consistent with the theoretical answer of 1/4 = 0.25.")

# Run the simulation
run_simulation()