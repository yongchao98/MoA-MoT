import math
import random

def solve_probability():
    """
    Calculates the probability that a circle of radius 6, thrown randomly onto a
    unit grid, intersects exactly 47 cells, using a Monte Carlo simulation.
    """
    R = 6.0
    R2 = R * R
    TARGET_DIFFERENCE = 22  # This corresponds to intersecting 47 cells.

    num_samples = 5_000_000
    success_count = 0

    # The problem asks for the number of cells intersected by the circumference.
    # Our analysis shows this leads to the condition N_xy - N_int = 22.
    # N_xy: number of grid vertices inside the circle.
    # N_int: number of grid cells completely inside the circle.

    for _ in range(num_samples):
        # Generate a random center (x0, y0) in the unit square [0, 1] x [0, 1]
        x0 = random.random()
        y0 = random.random()

        # The range of integer coordinates k, l to check for vertices and cells.
        # For a center in [0,1]x[0,1] and R=6, vertices can have k in [-5, 6], l in [-5, 6].
        # A safe check range is k, l in [-6, 7].
        k_min_check = -6
        k_max_check = 7
        l_min_check = -6
        l_max_check = 7

        n_xy = 0
        interior_vertices = set()
        for k in range(k_min_check, k_max_check + 1):
            for l in range(l_min_check, l_max_check + 1):
                if (k - x0)**2 + (l - y0)**2 < R2:
                    n_xy += 1
                    interior_vertices.add((k, l))

        n_int = 0
        # A cell [k,k+1]x[l,l+1] is interior if all four vertices are interior.
        # The range for the bottom-left corner (k,l) of the cell must be adjusted.
        for k in range(k_min_check, k_max_check):
             for l in range(l_min_check, l_max_check):
                if (
                    (k, l) in interior_vertices and
                    (k + 1, l) in interior_vertices and
                    (k, l + 1) in interior_vertices and
                    (k + 1, l + 1) in interior_vertices
                ):
                    n_int += 1

        if (n_xy - n_int) == TARGET_DIFFERENCE:
            success_count += 1
            
    probability = success_count / num_samples

    print(f"The calculation is based on the formula: N_intersected = 25 + N_vertices - N_cells.")
    print(f"We need N_intersected = 47, so we search for cases where N_vertices - N_cells = 22.")
    print(f"\nSimulation Results:")
    print(f"Total number of random throws: {num_samples}")
    print(f"Number of throws intersecting 47 cells: {success_count}")
    print(f"\nThe probability is calculated by the equation:")
    # "output each number in the final equation!"
    print(f"P = {success_count} / {num_samples}")
    print(f"Approximate probability: {probability:.4g}")

solve_probability()
<<<0.02107>>>