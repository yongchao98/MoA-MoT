import random
import math

def solve_circumference_problem():
    """
    Calculates the probability that a circumference of radius 6, thrown at random
    onto a plane with a square grid of mesh 1, intersects exactly 47 grid cells.
    """
    R = 6.0
    R2 = R * R
    NUM_SAMPLES = 10_000_000
    count_for_47_cells = 0

    # The range of cell indices to check.
    # For a center (xc, yc) in [0,1)x[0,1) and R=6, the circle can span from
    # xc-R (in [-6,-5)) to xc+R (in [6,7)).
    # So, integer indices i and j will range from -6 to 6.
    i_range = range(-6, 7)
    j_range = range(-6, 7)

    for _ in range(NUM_SAMPLES):
        # Generate a random center (xc, yc) in the unit square [0, 1) x [0, 1)
        xc = random.random()
        yc = random.random()

        intersected_cells = 0
        for i in i_range:
            for j in j_range:
                # For each cell [i, i+1] x [j, j+1], check for intersection.

                # Calculate the squared minimum distance from the center to the cell.
                # This is the distance to the closest point in the cell.
                dx_min = 0.0
                if xc < i:
                    dx_min = i - xc
                elif xc > i + 1:
                    dx_min = xc - (i + 1)

                dy_min = 0.0
                if yc < j:
                    dy_min = j - yc
                elif yc > j + 1:
                    dy_min = yc - (j + 1)

                d2_min = dx_min**2 + dy_min**2

                # If the minimum distance is greater than R, the circle is entirely
                # outside the cell.
                if d2_min > R2:
                    continue

                # Calculate the squared maximum distance from the center to the cell.
                # This is the distance to the furthest corner of the cell.
                dx_max = max(abs(i - xc), abs(i + 1 - xc))
                dy_max = max(abs(j - yc), abs(j + 1 - yc))
                d2_max = dx_max**2 + dy_max**2

                # If the maximum distance is less than R, the cell is entirely
                # inside the circle.
                if d2_max < R2:
                    continue

                # If we reach here, the cell is intersected by the circumference.
                intersected_cells += 1

        if intersected_cells == 47:
            count_for_47_cells += 1

    probability = count_for_47_cells / NUM_SAMPLES

    print(f"Radius (R): {R}")
    print(f"Total number of random samples: {NUM_SAMPLES}")
    print(f"Number of times the circumference intersected exactly 47 cells: {count_for_47_cells}")
    print(f"The probability is the ratio of these two numbers.")
    print(f"P(N=47) = {count_for_47_cells} / {NUM_SAMPLES}")
    print(f"Approximate solution: {probability:.4f}")

solve_circumference_problem()