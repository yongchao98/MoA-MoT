import math
import random

def solve_circumference_problem():
    """
    Calculates the probability that a circle of radius 6, thrown randomly onto a 
    1x1 grid, intersects exactly 47 cells, using a Monte Carlo simulation.
    """
    R = 6.0
    R2 = R * R
    N_TARGET = 47
    # A large number of trials is needed for an accurate approximation.
    N_TRIALS = 10_000_000

    count_target = 0

    # Define the range of grid cells to check for each trial.
    # A cell [n, n+1]x[m, m+1] can only be intersected if its minimum distance
    # to the circle's center is less than or equal to R.
    # Since the center (x_c, y_c) is in [0,1]x[0,1], the relevant cell indices
    # n and m range from -7 to 7.
    n_min, n_max = -7, 7
    m_min, m_max = -7, 7

    for i in range(N_TRIALS):
        # 1. Generate a random center (x_c, y_c) in the unit square [0, 1] x [0, 1]
        x_c = random.random()
        y_c = random.random()

        # 2. Count the number of intersected cells
        intersected_cells = 0
        for n in range(n_min, n_max + 1):
            for m in range(m_min, m_max + 1):
                # 3. Check the intersection condition for the cell [n, n+1] x [m, m+1]

                # Calculate the squared minimum distance (d_min^2) from the center 
                # (x_c, y_c) to the closest point in the cell.
                dx_min = 0
                if x_c < n:
                    dx_min = n - x_c
                elif x_c > n + 1:
                    dx_min = x_c - (n + 1)

                dy_min = 0
                if y_c < m:
                    dy_min = m - y_c
                elif y_c > m + 1:
                    dy_min = y_c - (m + 1)

                d2_min = dx_min**2 + dy_min**2

                # If the circle is entirely outside the cell, skip to the next cell.
                if d2_min > R2:
                    continue

                # Calculate the squared maximum distance (d_max^2) from the center 
                # (x_c, y_c) to the furthest point in the cell (one of its corners).
                dx_max = max(abs(n - x_c), abs(n + 1 - x_c))
                dy_max = max(abs(m - y_c), abs(m + 1 - y_c))
                d2_max = dx_max**2 + dy_max**2

                # If the cell is entirely contained within the circle, the circumference
                # does not intersect it. Skip to the next cell.
                if d2_max < R2:
                    continue

                # If d_min <= R <= d_max, the circumference intersects the cell.
                intersected_cells += 1

        # 4. Check if the number of intersected cells is the target number.
        if intersected_cells == N_TARGET:
            count_target += 1

    # 5. Calculate the final probability.
    probability = count_target / N_TRIALS

    print("--- Monte Carlo Simulation Results ---")
    print(f"Circle Radius (R): {R}")
    print(f"Target number of intersected cells: {N_TARGET}")
    print(f"Total trials performed: {N_TRIALS}")
    print(f"Number of successful trials (hits): {count_target}")
    print("\nThe probability is the ratio of hits to total trials:")
    print(f"P(N=47) = {count_target} / {N_TRIALS}")
    print(f"Approximate Probability = {probability:.4g}") # .4g for 4 significant digits

if __name__ == '__main__':
    solve_circumference_problem()