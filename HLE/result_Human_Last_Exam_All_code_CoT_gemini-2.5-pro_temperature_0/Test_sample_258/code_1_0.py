import math

def solve():
    """
    Calculates the minimal and maximal numbers of grid cells a circle can cross.
    """
    R = 500

    # Case 1: Maximal number of cells
    # This occurs when the center is at (k + 0.5, l + 0.5).
    # We can set k=0, l=0 for simplicity, so center is (0.5, 0.5).
    x0_max = 0.5
    y0_max = 0.5

    # Number of vertical lines crossed, |n - x0| < R
    # n is in (x0 - R, x0 + R) = (0.5 - 500, 0.5 + 500) = (-499.5, 500.5)
    # Integers n are from -499 to 500.
    n_min_max = -499
    n_max_max = 500
    N_v_lines_max = n_max_max - n_min_max + 1

    # Number of horizontal lines crossed is the same by symmetry.
    N_h_lines_max = N_v_lines_max

    # The number of cells crossed is the total number of intersections.
    # Each line is crossed twice.
    max_cells = 2 * N_v_lines_max + 2 * N_h_lines_max
    
    print("Calculation for Maximal number of cells:")
    print(f"Radius R = {R}")
    print(f"Center position (x0, y0) is chosen to be ({x0_max}, {y0_max}) to maximize crossings.")
    print(f"The vertical lines x=n are crossed for n in ({x0_max - R}, {x0_max + R}), which is ({x0_max - R}, {x0_max + R}).")
    print(f"The integer values for n range from {n_min_max} to {n_max_max}.")
    print(f"Number of vertical lines crossed (N_v_lines) = {n_max_max} - {n_min_max} + 1 = {N_v_lines_max}")
    print(f"Number of horizontal lines crossed (N_h_lines) = {N_h_lines_max}")
    print(f"Maximal number of cells = 2 * N_v_lines + 2 * N_h_lines = 2 * {N_v_lines_max} + 2 * {N_h_lines_max} = {max_cells}")
    print("-" * 20)

    # Case 2: Minimal number of cells
    # This occurs when the center is very close to a grid point, e.g., (eps, eps).
    # We can let eps be a very small number, effectively x0=0, y0=0 for floor/ceil calculations.
    x0_min_repr = 0.000001 # a small epsilon
    y0_min_repr = 0.000001

    # Number of vertical strips crossed
    # N_v_strips = floor(x0 + R) - floor(x0 - R)
    N_v_strips_min = math.floor(x0_min_repr + R) - math.floor(x0_min_repr - R)

    # Number of horizontal strips crossed
    N_h_strips_min = math.floor(y0_min_repr + R) - math.floor(y0_min_repr - R)

    # The number of cells is the perimeter of the bounding box of strips.
    min_cells = 2 * (N_v_strips_min + N_h_strips_min - 2)

    print("Calculation for Minimal number of cells:")
    print(f"Radius R = {R}")
    print(f"Center position (x0, y0) is chosen to be close to a grid point, e.g., ({x0_min_repr}, {y0_min_repr}).")
    print(f"Number of vertical strips crossed (N_v_strips) = floor({x0_min_repr + R}) - floor({x0_min_repr - R}) = {N_v_strips_min}")
    print(f"Number of horizontal strips crossed (N_h_strips) = floor({y0_min_repr + R}) - floor({y0_min_repr - R}) = {N_h_strips_min}")
    print(f"Minimal number of cells = 2 * (N_v_strips + N_h_strips - 2) = 2 * ({N_v_strips_min} + {N_h_strips_min} - 2) = {min_cells}")
    print("-" * 20)
    
    print(f"The final answer (min, max) is ({min_cells}, {max_cells}).")

solve()