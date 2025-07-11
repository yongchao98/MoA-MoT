import math
import random

def solve_circumference_probability():
    """
    Calculates the probability that a circumference of radius 6, thrown at random
    onto a plane with a square grid of mesh 1, intersects exactly 47 grid cells.
    This is done using a Monte Carlo simulation.
    """
    R = 6.0
    R_sq = R**2
    TARGET_INTERSECTIONS = 47
    # A large number of trials is needed for an accurate result.
    # This may take a few seconds to run.
    NUM_TRIALS = 2_000_000

    # Set a seed for reproducibility of the random results.
    random.seed(123)

    success_count = 0

    # Define the integer range of cells to check.
    # A circle with center in [0,1]x[0,1] and radius 6 can only intersect
    # cells with indices i in [-7, 6] and j in [-7, 6].
    i_min, i_max = -7, 7  # Loop from -7 to 6 inclusive
    j_min, j_max = -7, 7

    for _ in range(NUM_TRIALS):
        # 1. Generate a random center for the circle in the unit square.
        x_c = random.uniform(0, 1)
        y_c = random.uniform(0, 1)

        intersections = 0
        # 2. Iterate through all potentially intersected cells.
        for i in range(i_min, i_max):
            for j in range(j_min, j_max):
                # A cell is the square [i, i+1] x [j, j+1]

                # Find the squared distance from the center (x_c, y_c) to the
                # closest point in the cell.
                px = max(i, min(x_c, i + 1))
                py = max(j, min(y_c, j + 1))
                d_min_sq = (x_c - px)**2 + (y_c - py)**2

                # If the circle's radius is smaller than the minimum distance,
                # the circle is completely outside the cell.
                if d_min_sq >= R_sq:
                    continue

                # Find the squared distance to the farthest corner of the cell.
                d_max_sq = max(
                    (x_c - i)**2 + (y_c - j)**2,
                    (x_c - (i + 1))**2 + (y_c - j)**2,
                    (x_c - i)**2 + (y_c - (j + 1))**2,
                    (x_c - (i + 1))**2 + (y_c - (j + 1))**2
                )

                # If the radius is larger than the maximum distance, the cell
                # is completely inside the circle. The circumference doesn't
                # intersect its interior.
                # The condition d_min < R < d_max means the circumference
                # must pass through the cell.
                if d_max_sq > R_sq:
                    intersections += 1

        # 3. Check if the number of intersections matches the target.
        if intersections == TARGET_INTERSECTIONS:
            success_count += 1

    # 4. Calculate and print the final probability.
    probability = success_count / NUM_TRIALS

    print("--- Monte Carlo Simulation Results ---")
    print(f"Radius (R): {R}")
    print(f"Target Intersections: {TARGET_INTERSECTIONS}")
    print(f"Total Trials: {NUM_TRIALS}")
    print(f"Successful Trials (intersections == {TARGET_INTERSECTIONS}): {success_count}")
    print("\nFinal Probability Equation:")
    print(f"P(N=47) ≈ {success_count} / {NUM_TRIALS}")
    print(f"\nApproximate probability: {probability:.4g}")

# Run the simulation
solve_circumference_probability()