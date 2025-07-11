import random
import math

def solve():
    """
    Calculates the probability that the number of intersected cells is 47
    using a Monte Carlo simulation.
    """
    # Parameters for the simulation
    R = 6.0
    R2 = R * R
    # A large number of trials is needed for precision.
    num_trials = 5 * 10**6
    # Use a fixed seed for reproducibility of the result.
    random.seed(0)

    def get_min_dist_sq(i, j, xc, yc):
        """Calculates the square of the minimum distance from (xc, yc) to cell Q_{i,j}."""
        dx = 0.0
        if xc < i:
            dx = i - xc
        elif xc > i + 1:
            dx = xc - (i + 1)

        dy = 0.0
        if yc < j:
            dy = j - yc
        elif yc > j + 1:
            dy = yc - (j + 1)

        return dx*dx + dy*dy

    def get_max_dist_sq(i, j, xc, yc):
        """Calculates the square of the maximum distance from (xc, yc) to cell Q_{i,j}."""
        # The maximum distance is to one of the four corners of the cell
        dx1 = xc - i
        dx2 = xc - (i + 1)
        dy1 = yc - j
        dy2 = yc - (j + 1)
        return max(dx1*dx1 + dy1*dy1, dx1*dx1 + dy2*dy2, dx2*dx2 + dy1*dy1, dx2*dx2 + dy2*dy2)

    def count_intersections(xc, yc):
        """Counts the number of grid cells intersected by the circle."""
        count = 0
        # A cell Q_{i,j} can be intersected only if it's close to the circle.
        # The loop range covers all potentially intersected cells for R=6 and xc,yc in [0,1].
        # i must be in floor(xc-R-1) to ceil(xc+R) => [-7, 7]
        # j must be in floor(yc-R-1) to ceil(yc+R) => [-7, 7]
        for i in range(-7, 8):
            for j in range(-7, 8):
                # An intersection occurs if the cell is not fully inside or fully outside the circle.
                # This is checked by comparing min and max distances from the center to the cell with R.
                min_d2 = get_min_dist_sq(i, j, xc, yc)
                if min_d2 <= R2:
                    max_d2 = get_max_dist_sq(i, j, xc, yc)
                    if max_d2 >= R2:
                        count += 1
        return count

    # Main Monte Carlo loop
    count47 = 0
    for _ in range(num_trials):
        # Sample a random center (xc, yc) uniformly from the unit square [0, 1] x [0, 1]
        xc = random.random()
        yc = random.random()

        n_intersect = count_intersections(xc, yc)
        if n_intersect == 47:
            count47 += 1

    probability = count47 / num_trials

    # Output the components of the probability calculation
    print(f"Total trials: {num_trials}")
    print(f"Favorable outcomes (number of intersections is 47): {count47}")
    print(f"The final equation is: Probability = {count47} / {num_trials}")
    
    # Output the final approximate probability
    print(f"Approximate probability: {probability:.4g}")

solve()