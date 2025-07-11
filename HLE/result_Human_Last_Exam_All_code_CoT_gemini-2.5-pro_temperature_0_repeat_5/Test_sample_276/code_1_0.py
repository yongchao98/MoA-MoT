import random
import math

def solve_circumference_problem():
    """
    Calculates the probability that a circumference of radius 6, thrown at random
    onto a unit square grid, intersects exactly 47 cells.
    """
    
    # --- Plan ---
    # 1. Define parameters for the simulation.
    # 2. Create a function to count cell intersections for a given circle center.
    # 3. Run a Monte Carlo simulation by generating random circle centers.
    # 4. Count successful trials (where intersections == 47).
    # 5. Calculate and print the probability.
    
    # Step 1: Define parameters
    R = 6.0
    R_sq = R * R
    TARGET_INTERSECTIONS = 47
    NUM_SAMPLES = 2_000_000 # Using a large number for better accuracy

    def count_intersections(x, y):
        """
        Counts the number of grid cells intersected by a circle's circumference.
        The circle has radius R and its center (x, y) is in [0, 1] x [0, 1].
        """
        count = 0
        # The circle's x-extent is [x-R, x+R]. Since x is in [0,1], this is within [-6, 7].
        # The integer coordinates 'i' of the bottom-left of potentially intersected cells
        # range from floor(x-R) to floor(x+R). This is from -6 to 7.
        # So we check cells with i from -6 to 7, and j from -6 to 7.
        for i in range(-6, 8):
            for j in range(-6, 8):
                # The cell is the square [i, i+1] x [j, j+1]

                # Find the squared distance from the circle's center (x,y)
                # to the closest point in the cell.
                px = max(i, min(x, i + 1))
                py = max(j, min(y, j + 1))
                min_dist_sq = (x - px)**2 + (y - py)**2

                # If the closest point is further than R, the whole cell is outside the circle.
                if min_dist_sq > R_sq:
                    continue

                # Find the squared distance from the circle's center (x,y)
                # to the farthest point in the cell. This must be one of the 4 corners.
                dx_max = max(abs(i - x), abs(i + 1 - x))
                dy_max = max(abs(j - y), abs(j + 1 - y))
                max_dist_sq = dx_max**2 + dy_max**2

                # If the farthest point is closer than R, the whole cell is inside the circle.
                if max_dist_sq < R_sq:
                    continue
                
                # If we are here, the cell is neither completely inside nor completely outside.
                # Therefore, the circumference intersects it.
                count += 1
        return count

    # Step 3 & 4: Run Monte Carlo simulation
    success_count = 0
    for _ in range(NUM_SAMPLES):
        # Generate a random center (x, y) in the unit square [0, 1] x [0, 1]
        rand_x = random.random()
        rand_y = random.random()

        num_intersected = count_intersections(rand_x, rand_y)

        if num_intersected == TARGET_INTERSECTIONS:
            success_count += 1

    # Step 5: Calculate and print the probability
    probability = success_count / NUM_SAMPLES

    print(f"A circumference of radius R = {R} is considered.")
    print(f"The target number of intersected cells is {TARGET_INTERSECTIONS}.")
    print(f"The simulation was run with {NUM_SAMPLES} samples.")
    print(f"The number of times the circumference intersected exactly {TARGET_INTERSECTIONS} cells was: {success_count}.")
    print(f"The final probability equation is: P = {success_count} / {NUM_SAMPLES}")
    # Format to 4 significant digits
    print(f"The approximate probability is: {probability:.4g}")

solve_circumference_problem()
<<<0.0236>>>